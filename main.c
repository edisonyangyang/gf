//--------------------------------------
//@author:edison
//@email: edyyoung@163.com


//--------------------------
#include "stdio.h"

//gf求解方程组的row,col
//todo: Adjust this to a variable parameter 
#define RowM 8
#define ColN 8 

//-------------------------------补充函数--------------------------------------------------
//有限域计算，gfmath
unsigned char  gf_Logtable[256] =
{ 0, 0, 25, 1, 50, 2, 26, 198, 75, 199, 27, 104, 51, 238, 223, 3,
100, 4, 224, 14, 52, 141, 129, 239, 76, 113, 8, 200, 248, 105, 28, 193,
125, 194, 29, 181, 249, 185, 39, 106, 77, 228, 166, 114, 154, 201, 9, 120,
101, 47, 138, 5, 33, 15, 225, 36, 18, 240, 130, 69, 53, 147, 218, 142,
150, 143, 219, 189, 54, 208, 206, 148, 19, 92, 210, 241, 64, 70, 131, 56,
102, 221, 253, 48, 191, 6, 139, 98, 179, 37, 226, 152, 34, 136, 145, 16,
126, 110, 72, 195, 163, 182, 30, 66, 58, 107, 40, 84, 250, 133, 61, 186,
43, 121, 10, 21, 155, 159, 94, 202, 78, 212, 172, 229, 243, 115, 167, 87,
175, 88, 168, 80, 244, 234, 214, 116, 79, 174, 233, 213, 231, 230, 173, 232,
44, 215, 117, 122, 235, 22, 11, 245, 89, 203, 95, 176, 156, 169, 81, 160,
127, 12, 246, 111, 23, 196, 73, 236, 216, 67, 31, 45, 164, 118, 123, 183,
204, 187, 62, 90, 251, 96, 177, 134, 59, 82, 161, 108, 170, 85, 41, 157,
151, 178, 135, 144, 97, 190, 220, 252, 188, 149, 207, 205, 55, 63, 91, 209,
83, 57, 132, 60, 65, 162, 109, 71, 20, 42, 158, 93, 86, 242, 211, 171,
68, 17, 146, 217, 35, 32, 46, 137, 180, 124, 184, 38, 119, 153, 227, 165,
103, 74, 237, 222, 197, 49, 254, 24, 13, 99, 140, 128, 192, 247, 112, 7 };

unsigned char gf_Alogtable[256] =
{ 1, 3, 5, 15, 17, 51, 85, 255, 26, 46, 114, 150, 161, 248, 19, 53,
95, 225, 56, 72, 216, 115, 149, 164, 247, 2, 6, 10, 30, 34, 102, 170,
229, 52, 92, 228, 55, 89, 235, 38, 106, 190, 217, 112, 144, 171, 230, 49,
83, 245, 4, 12, 20, 60, 68, 204, 79, 209, 104, 184, 211, 110, 178, 205,
76, 212, 103, 169, 224, 59, 77, 215, 98, 166, 241, 8, 24, 40, 120, 136,
131, 158, 185, 208, 107, 189, 220, 127, 129, 152, 179, 206, 73, 219, 118, 154,
181, 196, 87, 249, 16, 48, 80, 240, 11, 29, 39, 105, 187, 214, 97, 163,
254, 25, 43, 125, 135, 146, 173, 236, 47, 113, 147, 174, 233, 32, 96, 160,
251, 22, 58, 78, 210, 109, 183, 194, 93, 231, 50, 86, 250, 21, 63, 65,
195, 94, 226, 61, 71, 201, 64, 192, 91, 237, 44, 116, 156, 191, 218, 117,
159, 186, 213, 100, 172, 239, 42, 126, 130, 157, 188, 223, 122, 142, 137, 128,
155, 182, 193, 88, 232, 35, 101, 175, 234, 37, 111, 177, 200, 67, 197, 84,
252, 31, 33, 99, 165, 244, 7, 9, 27, 45, 119, 153, 176, 203, 70, 202,
69, 207, 74, 222, 121, 139, 134, 145, 168, 227, 62, 66, 198, 81, 243, 14,
18, 54, 90, 238, 41, 123, 141, 140, 143, 138, 133, 148, 167, 242, 13, 23,
57, 75, 221, 124, 132, 151, 162, 253, 28, 36, 108, 180, 199, 82, 246, 1 };

//加法
unsigned char gf_add(unsigned char a, unsigned char b)
{
	return a ^ b;
}

//乘法：
unsigned char  gf_mul(unsigned char a, unsigned char b) 
{
	if (a && b)
		return gf_Alogtable[(gf_Logtable[a] + gf_Logtable[b]) % 255];
	else return 0;
}



//除法实现
unsigned char gf_div(unsigned char a, unsigned char b) 
{
	int j;
	if (a == 0)
		return (0);
	if ((j = gf_Logtable[a] - gf_Logtable[b]) < 0) 
		j += 255;
	return (gf_Alogtable[j]);
}
//逆元
unsigned char gf_inv(unsigned char in) 
{
	/* 0 is self inverting */
	if (in == 0)
		return 0;
	else
		return gf_Alogtable[(255 - gf_Logtable[in])];
}

// 功能：有限域矩阵相乘  
// 形参：(输入)矩阵A，矩阵A行数row和列数col，矩阵B，矩阵B行数row和列数col  
// 返回：A*B 
unsigned char gf_matMul(unsigned char* A, int Arow, int Acol, unsigned char* B, int Brow, int Bcol,unsigned char *Temp)
{
	if (Acol != Brow)
	{
		printf("矩阵不可乘!\n");
		return NULL;
	}
	if (Acol == Brow)
	{
		int i, j, k;
		for (i = 0; i < Arow; i++)
		{
			for (j = 0; j < Bcol; j++)
			{
				Temp[Bcol*i + j] = 0;
				for (k = 0; k < Acol; k++)
					Temp[Bcol*i + j] = gf_add(Temp[Bcol*i + j], gf_mul(A[Acol*i + k], B[Bcol*k + j]));
			}
		}
	}
	return 1;
}
//------------------------------高斯消元法,gf_guass---------------------------------------
unsigned char gf_guass_RowM_ColN(unsigned char coef[],int coef_row,int coef_col, unsigned char b[],unsigned char x[])
{
	int i, j, k;
	unsigned char mi;
	unsigned char mx, tmp;
	unsigned char a[RowM][ColN]={0};

	for (i = 0; i < coef_row; i++)
	{
		for (j = 0; j < coef_col; j++)
			a[i][j] = coef[i * coef_col + j];
	}

	for (i = 0; i<coef_col - 1; i++)
	{
		//找主元素 
		for (j = i + 1, mi = i, mx = a[i][i]; j<coef_row; j++)
		{
			if (a[j][i]>mx)
			{
				mi = j;
				mx = a[j][i];
			}
		}

		//交换两行  
		if (i<mi)
		{
			tmp = b[i];
			b[i] = b[mi];
			b[mi] = tmp;
			for (j = i; j< coef_col; j++)
			{
				tmp = a[i][j];
				a[i][j] = a[mi][j];
				a[mi][j] = tmp;
			}
		}

		//高斯消元 
		for (j = i + 1; j<coef_row; j++)
		{
			if(a[i][i]==0)
				return -1;
			else
			{
				tmp = gf_div(a[j][i], a[i][i]);
				b[j] = gf_add(b[j], gf_mul(b[i], tmp));
				for (k = i; k<coef_col; k++)
					a[j][k] = gf_add(a[j][k], gf_mul(a[i][k], tmp));
			}
		}
	}
	//求解方程
	//判断是否为0，即是否能偶求解
	if (a[coef_col - 1][coef_col - 1] == 0)
	{
		return 2;
	}
	else
	{
		x[coef_col - 1] = gf_div(b[coef_col - 1], a[coef_col - 1][coef_col - 1]);
		for (i = coef_col - 2; i >= 0; i--)
		{
			x[i] = b[i];
			for (j = i + 1; j <coef_col; j++)
				x[i] = gf_add(x[i], gf_mul(a[i][j], x[j]));
			x[i] = gf_div(x[i], a[i][i]);
		}
		return 1;	
	}

}










