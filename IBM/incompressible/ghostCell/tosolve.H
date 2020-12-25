#ifndef TOSOLVE_H
#define TOSOLVE_H

//// 变量篇
// Fluid cells indicator, marking only live fluid cells
// gammaPtr

// fluid cells indicator, marking live and IB cells
// gammaExtPtr

// List of fluid cells next to immersed boundary (IB cells)
// ibCellsPtr

// 根据单元格所在的面，寻找周围的单元格
// maybe ibFacePtr or ibFaceCellsPtr

// Normals at IB points

// 1. 标记流体网格，死网格和鬼单元格
// 2. 寻找距离鬼单元格最近的流体网格

//// 函数篇
// 1. void makeGamma() const


// 2. void makeIbcells() const

// 3. void makeDeadCells() const

// 4. void makeIbNormals() const

// 5. void make IbCellCells() const

// 6. label findNearestCell(const point& location) const;

// clear function

// 7. dirichlet边界条件设置：那个矩阵，求解系数，如何存储系数，读取速度场

// 8. norman边界条件设置：

// 9. 如何计算颗粒所受的力，根据速度场积分把。

#endif // TOSOLVE_H

#include <iostream>
#include <math.h>	// fabs
#include <stdio.h>	// pintf

using namespace std;

// 交换第i行和第j行
// 这里包括了b那一列
void swapA(int i, int j, int n, double A[][5]){
  for(int x = 0; x <= n; x++){
    double temp = A[i][x];
    A[i][x] = A[j][x];
    A[j][x] = temp;
  }
}

void getResult(int n, double e, int N, double A[][5], double re[])
{
  int i,j,k;
  // n个未知参数
  for(i = 0; i < n; i++){
    re[i] = 0.0;
  }

  // 先检验对角线元素是否为0，如果为0则交换某两列
  // 什么时候会交换失败？
  for(i = 0; i < n; i++){
    if(fabs(A[i][i])<=1e-2){
      for(j = 0; j < n; j++){
    if(fabs(A[j][i]) > 1e-1){
      swapA(i,j,n,A);
      break;
    }
      }
      if(j == n){
        printf("系数矩阵不合法! \n");
      }
      i = 0;	// 每次从头找
    }
  }
  for(i = 0; i < n; i++){
    for(j = 0; j < n+1; j++){
      printf("%lf ", A[i][j]);
    }
    cout << endl;
  }

  // iteration
  k = 0;
  double x[10];
  // initialize x
  for(i = 0; i < n; i++){
    x[i] = 0.0;
  }
  while(k <= N){
    k++;
    if(k > N){
      printf("迭代失败!\n");
      break;
    }
    for(i = 0; i < n; i++){
      re[i] = A[i][n];
      for(j = 0; j < n; j++){
    if(j!=i){
      re[i] = re[i] - A[i][j]*x[j];
    }
      }
      re[i] = re[i] / A[i][i];
    }
    // 当最大的x误差小于e则退出
    double maxXerror = 0.0;
    for(i = 0; i < n; i++){
      if(fabs(x[i]-re[i]) > maxXerror){
        maxXerror = fabs(x[i] - re[i]);
      }
    }
    if(maxXerror < e){
      return;
    }
    printf("第%d步迭代结果: ", k);
    for(i = 0; i < n; i++){
      printf("%lf ", re[i]);
    }
    cout << endl;
    // 否则，继续
    for(i = 0; i < n; i++){
      x[i] = re[i];
    }
  }
}

int main()
{
  int x,y;

    printf("-------jacobi iteration--------\n\n");
  cout << "输入未知数个数与方程个数: ";
  cin >> x >> y;
  if(x!=y){
    cout << "超定或欠定! " << endl;
    return 0;
  }

  double A[4][5];
  double re[4];

  // 输入增广矩阵
  printf("输入增广矩阵:\n");
  for(int i = 0; i < x; i++){
    for(int j = 0; j < x+1; j++){
      cin >> A[i][j];
    }
  }
  cout << "输入精度e和最大迭代数: ";
  double e;
  int N;
  cin >> e >> N;
  getResult(x,e,N,A,re);
  for(int j = 0; j < x; j++){
    cout << re[j] << " ";
  }
  return 0;
}

// 1. 这个fluidPtr需不需要在求解器中声明

