#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {

    /**
    TODO:

    Complete the initialization. See ukf.h for other member properties.

    Hint: one or more values initialized above might be wildly off...
    */

  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initialize the state dimension
  n_x_ = 5;

  // initialize the augmented state dimension
  n_aug_ = 7;

  // initialize the sigma point spreading parameter
  lambda_ = 3 - n_aug_;

  //initialize number of times the Radar NIS exceeded the 95 % number
  radar_nis_exceed_num_ = 0;

  //initialize number of times the Lidar NIS exceeded the 95 % number
  lidar_nis_exceed_num_ = 0;

  // initial state vector
  x_ = VectorXd(n_x_);

  // initial covariance matrix
  P_ = MatrixXd(n_x_, n_x_);

  // Process noise standard deviation longitudinal acceleration in m/s^2.
  //This value is chosen to minimize the NIS exceed number and the RMSE. here we assume the maximum acceleration is 1 m/s^s
  std_a_ = 0.5;

  // Process noise standard deviation yaw acceleration in rad/s^2
  //This value is chosen to minimize the NIS exceed number and the RMSE.
  std_yawdd_ = M_PI / 4.0;

  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;


  // this flag shall be set to true after the first call of ProcessMeasurement()
  is_initialized_ = false;

  // initialize the predicted sigma points matrix
  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ +1);

  // initialize time when the state is true, in us
  time_us_ = 0;

  // initialize Weights of sigma points
  weights_ = VectorXd(2 * n_aug_ + 1);

}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Make sure you switch between lidar and radar
  measurements.
  */

    /*****************************************************************************
      *  Initialization
      ****************************************************************************/
       if (!is_initialized_) {
           /**
       TODO:
         * Initialize the state ekf_.x_ with the first measurement.
         * Create the covariance matrix.
         * Remember: you'll need to convert radar from polar to cartesian coordinates.
       */
           // first measurement
           cout << "UKF: " << endl;

           // initialize the state vector valuse
           x_(2) = 5.0;
           x_(3) = 0.0;
           x_(4) = 0.005;

           //Initialize the state covariance
           P_ = MatrixXd::Identity(n_x_, n_x_);

           P_(0,0) = 0.2;
           P_(1,1) = 0.2;
           P_(2,2) = 0.8;
           P_(3,3) = 0.5;
           P_(4,4) = 0.2;

           if (meas_package.sensor_type_ == MeasurementPackage::RADAR)
           {
               //Convert radar from polar to cartesian coordinates and initialize state.

               float rho = meas_package.raw_measurements_(0);
               float phi = meas_package.raw_measurements_(1);

               x_(0) = rho * cos(phi);
               x_(1) = rho * sin(phi);

           }
           else if (meas_package.sensor_type_ == MeasurementPackage::LASER)
           {

               x_(0) = meas_package.raw_measurements_(0);
               x_(1) = meas_package.raw_measurements_(1);
           }

           time_us_ = meas_package.timestamp_;

           // done initializing, no need to predict or update
           is_initialized_ = true;

           return;

       }

       /*****************************************************************************
        *  Prediction
        ****************************************************************************/


       //compute the time elapsed between the current and previous measurements
       float dt = (meas_package.timestamp_ - time_us_) / 1000000.0;	//dt - expressed in seconds
       time_us_ = meas_package.timestamp_;


       Prediction(dt);


       /*****************************************************************************
        *  Update
        ****************************************************************************/


       if (meas_package.sensor_type_ == MeasurementPackage::RADAR)
       {
           // Radar updates

           UpdateRadar(meas_package);

       }
       else
       {
           // Laser updates

           UpdateLidar(meas_package);
       }

       // print the output
       cout << "x_ = " << x_ << endl;
       cout << "P_ = " << P_ << endl;

       cout << "radar_nis_exceed"<< radar_nis_exceed_num_<<endl;
       cout << "lidar_nis_exceed"<< lidar_nis_exceed_num_<<endl;

}

/**
   * Calculates the weights vector and sigma points matrix
   * from the given state vector, state covariance matrix and lambda
   * @param state The state vector
   * @param cov The state covariance matrix
   * @param lambda Sigma point spreading parameter
   * @param weights The weights vector that shall be calculated
   * @param sigma_pts The sigma points matrix that shall be calculated
   */
void UKF::GetWeightsAndSigmaPts(const Eigen::VectorXd &state, const Eigen::MatrixXd &cov, double lambda,
                           Eigen::VectorXd &weights, Eigen::MatrixXd &sigma_pts)
{
    long n_x = state.size();

    //create square root matrix
    MatrixXd L = cov.llt().matrixL();

    //Calculate sigma points and the weights vector
    sigma_pts.col(0)  = state;
    weights(0) = lambda/(lambda+n_x);

    for (int i = 0; i< n_x; i++)
    {
      sigma_pts.col(i+1)       = state + sqrt(lambda + n_x) * L.col(i);
      sigma_pts.col(i+1+n_x) = state - sqrt(lambda + n_x) * L.col(i);

      weights(i+1) = 0.5/(n_x+lambda);
      weights(i+1+n_x) = 0.5/(n_x+lambda);
    }

}


/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
    /**
  TODO:

  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */

    //create augmented mean vector
    VectorXd x_aug = VectorXd(n_aug_);

    //create augmented state covariance
    MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);

    //create sigma point matrix
    MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);


    //create augmented mean state
    x_aug.head(n_x_) = x_;
    x_aug(5) = 0;
    x_aug(6) = 0;

    //create augmented covariance matrix
    P_aug.fill(0.0);
    P_aug.topLeftCorner(n_x_,n_x_) = P_;
    P_aug(5,5) = std_a_ * std_a_;
    P_aug(6,6) = std_yawdd_ * std_yawdd_;

    // Calculate the weights and sigma points from the state and state covariance
    GetWeightsAndSigmaPts (x_aug, P_aug, lambda_, weights_, Xsig_aug);

    //loop over the augmented simga points and calculate the predicted sigma points
     for (int i = 0; i< 2*n_aug_+1; i++)
     {
       //extract values for better readability
       double p_x = Xsig_aug(0,i);
       double p_y = Xsig_aug(1,i);
       double v = Xsig_aug(2,i);
       double psi = Xsig_aug(3,i);
       double psid = Xsig_aug(4,i);
       double nu_a = Xsig_aug(5,i);
       double nu_psidd = Xsig_aug(6,i);

       //declare the predicted sigma point elements
       double px_pred, py_pred, v_pred, psi_pred, psid_pred;

       //check to prevent division by zero
       if (fabs(psid) > 0.001) {
           px_pred = p_x + v/psid * ( sin (psi + psid*delta_t) - sin(psi));
           py_pred = p_y + v/psid * ( cos(psi) - cos(psi+psid*delta_t) );
       }
       else {
           px_pred = p_x + v*delta_t*cos(psi);
           py_pred = p_y + v*delta_t*sin(psi);
       }

       v_pred = v;
       psi_pred = psi + psid*delta_t;
       psid_pred = psid;

       //addition of the process noise
       px_pred = px_pred + 0.5*nu_a*delta_t*delta_t * cos(psi);
       py_pred = py_pred + 0.5*nu_a*delta_t*delta_t * sin(psi);
       v_pred = v_pred + nu_a*delta_t;
       psi_pred = psi_pred + 0.5*nu_psidd*delta_t*delta_t;
       psid_pred = psid_pred + nu_psidd*delta_t;

       // predicted sigma point
       Xsig_pred_(0,i) = px_pred;
       Xsig_pred_(1,i) = py_pred;
       Xsig_pred_(2,i) = v_pred;
       Xsig_pred_(3,i) = psi_pred;
       Xsig_pred_(4,i) = psid_pred;
     }

     // Calculate the predicted state mean
     x_.fill(0.0);
     for (int i = 0; i < 2 * n_aug_ + 1; i++) {
         x_ = x_+ weights_(i) * Xsig_pred_.col(i);
     }

     // Calculate the predicted state covariance matrix
     P_.fill(0.0);
     for (int i = 0; i < 2 * n_aug_ + 1; i++) {

         // state difference
         VectorXd x_diff = Xsig_pred_.col(i) - x_;

         //angle normalization
         x_diff(3) = atan2(sin(x_diff(3)), cos(x_diff(3)));

         P_ = P_ + weights_(i) * x_diff * x_diff.transpose() ;
     }

}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
    /* The Lidar measurement function is linear so UKF is not needed here. I shall use the normal Kalman Filter */

    VectorXd y;
    MatrixXd Pzz;

    MatrixXd H = MatrixXd(2, 5);
    H << 1, 0, 0, 0, 0,
         0, 1, 0, 0, 0;

    //measurement covariance matrix - laser
    MatrixXd R_laser = MatrixXd(2,2);

    R_laser << std_laspx_*std_laspx_, 0,
               0, std_laspy_*std_laspy_;

    //Calculate the innovation estimate
    y = meas_package.raw_measurements_ - H * x_;

    MatrixXd PHt = P_ * H.transpose();

    //Calculate the innovation covariance
    Pzz = H * PHt + R_laser;

    //Calculate Kalman Gain
    MatrixXd K = PHt * Pzz.inverse();

    //Calculate the fused state expectation
    x_ = x_ + K * y;

    //Calculate the fused state covariance
    long x_size = x_.size();
    MatrixXd I = MatrixXd::Identity(x_size, x_size);

    MatrixXd temp = (I - K * H)*P_;

    P_ = temp;


    double NIS_RADAR = y.transpose() * Pzz.inverse() * y;


    if (NIS_RADAR > 5.991)
    {
        lidar_nis_exceed_num_ += 1;
    }
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {

    //set measurement dimension, radar can measure r, phi, and r_dot
    long n_z = meas_package.raw_measurements_.size();

    //create matrix for sigma points in measurement space
    MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);

    //transform predicted sigma points into measurement space
    for (int i = 0; i < 2 * n_aug_ + 1; i++) {

        // extract values for better readibility
        double p_x = Xsig_pred_(0,i);
        double p_y = Xsig_pred_(1,i);
        double v  = Xsig_pred_(2,i);
        double psi = Xsig_pred_(3,i);

        double vx = v * cos(psi);
        double vy = v * sin(psi);

        // measurement model
        double rho = sqrt(p_x*p_x + p_y*p_y);
        if(fabs(rho) < 0.0001)
        {

            cout << "UpdateRadar () - Error - Division by Zero" << endl;
            Zsig(0,i) = rho;                                            //rho
            Zsig(1,i) = atan2(p_y,p_x);                                 //phi
            Zsig(2,i) = 0;                                              //rho_dot
        }
        else
        {
            Zsig(0,i) = rho;                                            //rho
            Zsig(1,i) = atan2(p_y,p_x);                                 //phi
            Zsig(2,i) = (p_x*vx + p_y*vy ) / rho;                       //rho_dot
        }

    }

    //mean predicted measurement
    VectorXd z_pred = VectorXd(n_z);
    z_pred.fill(0.0);
    for (int i=0; i < 2*n_aug_+1; i++) {
        z_pred = z_pred + weights_(i) * Zsig.col(i);
    }

    //measurement covariance matrix Pzz
    MatrixXd Pzz = MatrixXd(n_z,n_z);

    Pzz.fill(0.0);
    for (int i = 0; i < 2 * n_aug_ + 1; i++) {
        //residual
        VectorXd z_diff = Zsig.col(i) - z_pred;

        //angle normalization
        z_diff(1) = atan2(sin(z_diff(1)), cos(z_diff(1)));

        Pzz = Pzz + weights_(i) * z_diff * z_diff.transpose();
    }

    //add measurement noise covariance matrix
    MatrixXd R = MatrixXd(n_z,n_z);
    R <<    std_radr_*std_radr_, 0, 0,
            0, std_radphi_*std_radphi_, 0,
            0, 0,std_radrd_*std_radrd_;
    Pzz = Pzz + R;


    //create matrix for cross correlation Pxz
    MatrixXd Pxz = MatrixXd(n_x_, n_z);

    //calculate cross correlation matrix
    Pxz.fill(0.0);
    for (int i = 0; i < 2 * n_aug_ + 1; i++) {

        //residual
        VectorXd z_diff = Zsig.col(i) - z_pred;

        //angle normalization
        z_diff(1) = atan2(sin(z_diff(1)), cos(z_diff(1)));

        // state difference
        VectorXd x_diff = Xsig_pred_.col(i) - x_;

        //angle normalization
        x_diff(3) = atan2(sin(x_diff(3)), cos(x_diff(3)));

        Pxz = Pxz + weights_(i) * x_diff * z_diff.transpose();
    }

    //Kalman gain K;
    MatrixXd K = Pxz * Pzz.inverse();

    //residual
    VectorXd z_diff = meas_package.raw_measurements_ - z_pred;

    //angle normalization
    z_diff(1) = atan2(sin(z_diff(1)), cos(z_diff(1)));

    //update state mean and covariance matrix
    x_ = x_ + K * z_diff;
    P_ = P_ - K * Pzz * K.transpose();


    double NIS_RADAR = z_diff.transpose() * Pzz.inverse() * z_diff;


    if (NIS_RADAR > 7.815)
    {
        radar_nis_exceed_num_ += 1;
    }



}
