#include <stdio.h>
#include <stdlib.h>
#include <math.h> 
#include <string.h>
#include <ctime>

#define CombinLen 2      //SNP 组合个数
#define MaxSNPState 2    //SNP = 0,1 ... MaxSNPState
#define T_Value 0        //区分高感和低感的阈值
#define N_P     1000     //测试P值的循环次数

const char path[100] = "D:\\Documents\\Visual Studio 2013\\Projects\\MDR\\data\\MDR-SampleData.txt";
const char k_path[100] = "D:\\Documents\\Visual Studio 2013\\Projects\\MDR\\data\\K.txt";

struct SNPState
{
	int SNP1;
	int SNP2;
	int State;              //0:空缺 1:高感 2:低感
};

struct SNPCombin
{
	float Accurancy;
	int index_i;
	int index_j;
	SNPState State[(MaxSNPState + 1)*(MaxSNPState + 1)]; // 9
};

void main()
{
	const int data_Line = 400;
	const int data_Col = 21;
	const int K_num = 10; 
	//const int N_P = 1000;
	void KFoldCrossValid(int **data, const int row, const int col, SNPCombin *SNPData);
	void ExamSNP(int **data, const int row, const int col, SNPCombin SNPData, float *ExamAccuracy);
	void SNPTrainAndTest(int SNP1, int SNP2, int **data, const int row, const int col, const int k_num, float *Accuray);

	/*先略过确定文件行数和列数的判断，按照文件的大小赋值*/
	char data_Name[data_Col][10];
	int data[data_Line][data_Col];

	// 1. 读取数据
	FILE *fp;
	fp = fopen(path, "r");

	for (int i = 0; i < data_Col; i++)    //读取标题行
		fscanf(fp, "%s", data_Name[i]);

	for (int i = 0; i < data_Line; i++)   //读取数据
	{
		for (int j = 0; j < data_Col; j++)
		{
			fscanf(fp, "%d", &data[i][j]);
		}
	}
	fclose(fp);

	// 2. 训练数据
	//随机划分K个子集
	int data_index[data_Line];
	for (int i = 0; i < data_Line; i++)
	{
		data_index[i] = rand() % data_Line;
		for (int j = 0; j < i; j++)
		{
			while (data_index[i] == data_index[j])
			{
				data_index[i] = rand() % data_Line;
				j = 0;
			}
		}
		//printf("%d\n", data_index[i]);
	}

	FILE *fp2;
	fp2 = fopen(k_path, "w");
	int K[data_Line][data_Col];                       //K是划分子集后的数据，顺序存储了各个子集 
	for (int i = 0; i < data_Line; i++)
	{
		for (int j = 0; j < data_Col; j++)
		{
			K[i][j] = data[data_index[i]][j];
			fprintf(fp2, "%d\t", K[i][j]);
		}
		fprintf(fp2, "\n");
	}
	fclose(fp2);

	clock_t start, end;
	start = clock();

	//一个子集做测试集，其它K-1个子集做训练集。从训练集中挑选训练准确度最高的N个自变量组合
	//void KFoldCrossValid(int **data, const int row, const int col, SNPCombin *SNPData);
	//输入：训练集数据，数据行数，数据列数
	//输出：自变量组合  
	const int K_line = data_Line - data_Line / K_num;
	SNPCombin *SNPData_tmp, SNPData[K_num];
	SNPData_tmp = (struct SNPCombin*)malloc(sizeof(struct SNPCombin));              //结构体指针初始化
	for (int k = 0; k < K_num; k++)
	{
		//数据准备
		int K_data[K_line][data_Col];
		int tmp = 0;
		for (int i = 0; i < data_Line; i++)
		{
			if (i < k * data_Line / K_num || i >= (k + 1) * data_Line / K_num)
			{
				for (int j = 0; j < data_Col; j++)
				{
					K_data[tmp][j] = K[i][j];
				}
				tmp++;
			}
		}

		KFoldCrossValid((int **)K_data, K_line, data_Col, SNPData_tmp);
		SNPData[k] = *SNPData_tmp;
	}
	free(SNPData_tmp);

	//计算交叉一致性CVC
	//计算自变量组合个数num
	//int num;
	//long int tmp1, tmp2 = 1;
	//tmp1 = 1;
	//for (int i = 0; i < CombinLen; i++)
	//{
	//tmp1 *=  data_Col - 1 - i;
	//tmp2 *=  CombinLen - i;
	//}
	//num = tmp1 / CombinLen; 

	//计算各自变量组合被选择的次数
	int *SNP_Count;
	SNP_Count = (int *)malloc((data_Col - 1) * (data_Col - 1) * sizeof(int));
	memset(SNP_Count, 0, (data_Col - 1) * (data_Col - 1) * sizeof(int));
	for (int i = 0; i < K_num; i++)
	{
		int tmp;
		tmp = SNPData[i].index_i * (data_Col - 1) + SNPData[i].index_j;
		SNP_Count[tmp]++;
	}
	//选择符合条件的,且CVC最大的自变量组合
	int CVC_Threshold = K_num / 2;
	int max_SNP = 0;
	for (int i = 0; i < (data_Col - 1) * (data_Col - 1); i++)
	{
		if (SNP_Count[i] > SNP_Count[max_SNP])
		{
			max_SNP = i;
		}
	}
	int max_SNP1 = max_SNP / data_Col;
	int max_SNP2 = max_SNP % data_Col;
	if (SNP_Count[max_SNP] < CVC_Threshold)             //选择最多次的SNP组合的次数应该大于CVC阈值
	{
		printf("没有符合条件的自变量组合。\n");
		return;
	}
	else
		printf("自变量组合为：SNP%d  SNP%d\n", max_SNP1 + 1, max_SNP2 + 1);

	free(SNP_Count);

	//用测试子集来测试选定的自变量组合的准确率
	//void ExamSNP(int **data, const int row, const int col, SNPCombin SNPData, float ExamAccuracy);
	//输入：测试子集，测试子集的行数和列数，测试SNP组合的结构体
	//输出：测试准确率
	const int Exam_line = data_Line / K_num;
	int Exam_data[Exam_line][data_Col];
	float ExamAccuracy[K_num] = { 0 };              //存储各组的测试准确率。不是选择的组合不计算，准确率为0
	int exam_count = 0;
	for (int k = 0; k < K_num; k++)
	{
		int tmp = 0;
		if (SNPData[k].index_i == max_SNP1 && SNPData[k].index_j == max_SNP2)
		{
			// 数据准备
			for (int i = 0; i < data_Line; i++)
			{
				if (k * data_Line / K_num <= i && i < (k + 1) * data_Line / K_num)
				{
					for (int j = 0; j < data_Col; j++)
					{
						Exam_data[tmp][j] = K[i][j];
					}
					tmp++;
				}
			}

			ExamSNP((int **)Exam_data, data_Line / K_num, data_Col, SNPData[k], &ExamAccuracy[k]);
			exam_count++;
		}
	}

	float max_Accuracy = 0;
	float avg_Accuracy = 0;
	int max_index = 0;
	for (int k = 0; k < K_num; k++)
	{
		if (max_Accuracy < ExamAccuracy[k])
		{
			max_Accuracy = ExamAccuracy[k];
			max_index = k;
		}
		avg_Accuracy += ExamAccuracy[k];
	}
	avg_Accuracy = avg_Accuracy / exam_count;

	//计算P值
	int P_count = 0;
	float P;
	for (int n = 0; n < N_P; n++)
	{
		for (int i = 0; i < data_Line; i++)
		{
			int tmp = rand() % 2;            //对分值随机置换
			data[i][data_Col - 1] = tmp;
		}

		//对指定的自变量组合进行训练和测试
		//输入：自变量组合的索引，数据，行数，列数
		//输出：测试准确率
		//void SNPTrainAndTest(int SNP1, int SNP2, (int **)data, const int row, const int col,  const int k_num,float *Accuray);
		float Accuracy;
		SNPTrainAndTest(max_SNP1, max_SNP2, (int **)data, data_Line, data_Col, K_num, &Accuracy);
		if (Accuracy >= max_Accuracy)
		{
			P_count++;
		}
	}
	P = (float) P_count / N_P;

	end = clock();
	//double time = (double)(end - start) / CLOCKS_PER_SEC;
	int time = (double)(end - start);

	printf("测试最大准确率为：%4.2f%% 平均准确率为：%4.2f%%\n", max_Accuracy * 100, avg_Accuracy * 100);
	printf("\nIf-Then Rules:\n");
	for (int i = 0; i < (MaxSNPState + 1)*(MaxSNPState + 1); i++)
	{
		printf("IF SNP1 = %d AND SNP2 = %d THEN CLASSIFY AS %d\n",
			SNPData[max_index].State[i].SNP1,
			SNPData[max_index].State[i].SNP2,
			SNPData[max_index].State[i].State - 1);
	}

	printf("\nP值:   %f\n", P);
	printf("\nCPU计算耗时: %dms\n", time);
	printf("\n");
	return;
}

void KFoldCrossValid(int **data, const int row, const int col, SNPCombin *SNPData)
{
	//for (int i = 0; i < col; i++)
	//{
	//	printf("%d\t", *((int*)data + i));          //手动解析二维数组才能识别: *((int*)data + i*col + j) == data[i][j]
	//}



	//计算训练准确率
	float TrainAccuracy[200][3] = { 0 };                 //训练准确率，暂时用这个方法定义，个数应该为num = 190 个位点，
	SNPState State[200][(MaxSNPState + 1)*(MaxSNPState + 1)];  //存储分类矩阵
	int tmp1 = 0;                                             //tmp1:第几个自变量变量组合的index
	for (int i = 0; i < col - 2; i++)
	{
		for (int j = i + 1; j < col - 1; j++)
		{
			float Case1Sum[(MaxSNPState + 1)*(MaxSNPState + 1)] = { 0 };
			float Case0Sum[(MaxSNPState + 1)*(MaxSNPState + 1)] = { 0 };                      // =9 存储病历分值和、对照分值和

			int tmp2 = 0;                                    //tmp2:一个自变量组合的第几种形式的index
			//计算每一种SNP自变量组合的训练准确率
			for (int k = 0; k < row; k++)
			{
				//计算各组合的病历分值和及对照分值和
				tmp2 = *((int*)data + k*col + i) * (MaxSNPState + 1) + *((int*)data + k*col + j);
				if (*((int*)data + k*col + col - 1) == 1)
				{
					//Case1Sum[tmp2] += *((int*)data + k*col + col);       //最后一列数据
					Case1Sum[tmp2] ++;
				}
				else
				{
					//Case0Sum[tmp2] += *((int*)data + k*col + col);       //最后一列数据
					Case0Sum[tmp2] ++;
				}
			}

			//区分高感和低感，这里暂时考虑用差值大于0或者小于0来做   case0:病历 case1:对照
			float Accuracy = 0;
			float AFSum = 0;
			for (int m = 0; m < (MaxSNPState + 1)*(MaxSNPState + 1); m++)
			{
				int SNP1 = m / (MaxSNPState + 1);
				int SNP2 = m % (MaxSNPState + 1);
				if (Case0Sum[m] == 0 && Case1Sum[m] == 0)             //空缺
				{
					State[tmp1][m].SNP1 = SNP1;
					State[tmp1][m].SNP2 = SNP2;
					State[tmp1][m].State = 0;
				}
				else
				{
					if (Case0Sum[m] - Case1Sum[m] >= T_Value)         //高感
					{
						State[tmp1][m].SNP1 = SNP1;
						State[tmp1][m].SNP2 = SNP2;
						State[tmp1][m].State = 1;
						Accuracy += Case0Sum[m];
					}
					else                                              //低感
					{
						State[tmp1][m].SNP1 = SNP1;
						State[tmp1][m].SNP2 = SNP2;
						State[tmp1][m].State = 2;
						Accuracy += Case1Sum[m];
					}
				}
				AFSum += Case0Sum[m] + Case1Sum[m];
			}

			//根据病历分值和及计算分值和，计算训练准确率
			TrainAccuracy[tmp1][0] = Accuracy / AFSum;
			TrainAccuracy[tmp1][1] = i;         //存储位置，可以推算出自变量组合  
			TrainAccuracy[tmp1][2] = j;
			tmp1++;
		}
	}

	//选择准确率最大的自变量组合
	float MaxAccuracy = 0;
	int Max;
	for (int i = 0; i < tmp1; i++)
	{
		if (TrainAccuracy[i][0] > MaxAccuracy)
		{
			MaxAccuracy = TrainAccuracy[i][0];
			Max = i;
		}
	}
	//J[MaxIndexi * (MaxSNPState + 1) + MaxIndexj] ++;          //记录自变量组合被选中的次数
	SNPData->Accurancy = MaxAccuracy;
	SNPData->index_i = TrainAccuracy[Max][1];
	SNPData->index_j = TrainAccuracy[Max][2];
	for (int i = 0; i < (MaxSNPState + 1)*(MaxSNPState + 1); i++)
		SNPData->State[i] = State[Max][i];

	return;
}


void ExamSNP(int **data, const int row, const int col, SNPCombin SNPData, float *ExamAccuracy)
{
	float AccuracySum = 0;
	float TotalSum = 0;

	for (int i = 0; i < row; i++)
	{
		int SNP1 = *((int*)data + i*col + SNPData.index_i);
		int SNP2 = *((int*)data + i*col + SNPData.index_j);
		for (int j = 0; j < (MaxSNPState + 1)*(MaxSNPState + 1); j++)
		{
			if (SNP1 == SNPData.State[j].SNP1 && SNP2 == SNPData.State[j].SNP2)
			{
				if (*((int*)data + i*col + col) == SNPData.State[j].State - 1)          //病历(0)且为高感(1),对照(1)且为低感(2)
				{
					AccuracySum++;
					TotalSum++;
				}
				else
				{
					TotalSum++;
				}
			}
		}
	}

	*ExamAccuracy = AccuracySum / TotalSum;

	return;
}


void SNPTrainAndTest(int m_SNP1, int m_SNP2, int **data, const int row, const int col, const int k_num, float *m_Accuray)
{
	float Case1Sum[(MaxSNPState + 1)*(MaxSNPState + 1)] = { 0 };
	float Case0Sum[(MaxSNPState + 1)*(MaxSNPState + 1)] = { 0 };                      // =9 存储病历分值和、对照分值和

	for (int k = 0; k < row - row / k_num; k++)
	{
		//计算各组合的病历分值和及对照分值和
		int tmp2 = *((int*)data + k*col + m_SNP1) * (MaxSNPState + 1) + *((int*)data + k*col + m_SNP2);
		if (*((int*)data + k*col + col - 1) == 1)
		{
			Case1Sum[tmp2] ++;
		}
		else
		{
			Case0Sum[tmp2] ++;
		}
	}
	//区分高感和低感，这里暂时考虑用差值大于0或者小于0来做   case0:病历 case1:对照
	//float Accuracy = 0;
	//float AFSum = 0;
	SNPState State[(MaxSNPState + 1)*(MaxSNPState + 1)];  //存储分类矩阵
	for (int m = 0; m < (MaxSNPState + 1)*(MaxSNPState + 1); m++)
	{
		int SNP1 = m / (MaxSNPState + 1);
		int SNP2 = m % (MaxSNPState + 1);
		if (Case0Sum[m] == 0 && Case1Sum[m] == 0)             //空缺
		{
			State[m].SNP1 = SNP1;
			State[m].SNP2 = SNP2;
			State[m].State = 0;
		}
		else
		{
			if (Case0Sum[m] - Case1Sum[m] >= T_Value)         //高感
			{
				State[m].SNP1 = SNP1;
				State[m].SNP2 = SNP2;
				State[m].State = 1;
			}
			else                                              //低感
			{
				State[m].SNP1 = SNP1;
				State[m].SNP2 = SNP2;
				State[m].State = 2;
			}
		}
	}

	//测试
	float AccuracySum = 0;
	float TotalSum = 0;

	for (int i = row - row / k_num; i < row; i++)
	{
		int SNP1 = *((int*)data + i*col + m_SNP1);
		int SNP2 = *((int*)data + i*col + m_SNP2);
		for (int j = 0; j < (MaxSNPState + 1)*(MaxSNPState + 1); j++)
		{
			if (SNP1 == State[j].SNP1 && SNP2 == State[j].SNP2)
			{
				if (*((int*)data + i*col + col) == State[j].State - 1)          //病历(0)且为高感(1),对照(1)且为低感(2)
				{
					AccuracySum++;
					TotalSum++;
				}
				else
				{
					TotalSum++;
				}
			}
		}
	}
		
	*m_Accuray = AccuracySum / TotalSum;


	return;
}