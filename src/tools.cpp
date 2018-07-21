#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
  TODO:
    * Calculate the RMSE here.
  */
  
	VectorXd RMSE = VectorXd(4);
	RMSE << 0, 0, 0, 0;
	if (estimations.size() == 0)
	{
		cout << "estimations vector size : 0" << endl;
		return RMSE;
	}
	if (ground_truth.size() != estimations.size())
	{
		cout << "estimations vector size doesn't match ground truth vector size" << endl;
		return RMSE;
	}
	for (int i = 0; i < int(estimations.size()); ++i)
	{
		VectorXd residual = estimations[i] - ground_truth[i];
		VectorXd residual_2 = residual.array()*residual.array();
		RMSE = RMSE + residual_2;
	}
	RMSE = RMSE / estimations.size();
	RMSE = RMSE.array().sqrt();

	return RMSE;
}