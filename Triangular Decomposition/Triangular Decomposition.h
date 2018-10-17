#pragma once
#include <Eigen/Dense>

void DoolittleLU(Eigen::MatrixXd A, Eigen::MatrixXd& L, Eigen::MatrixXd& U);
void MainElementDoolittleLU(Eigen::MatrixXd A,Eigen::MatrixXd& B, Eigen::MatrixXd& L, Eigen::MatrixXd& U);
void CroutLU(Eigen::MatrixXd A, Eigen::MatrixXd& L, Eigen::MatrixXd& U);