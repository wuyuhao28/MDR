#include <stdio.h>
#include <stdlib.h>
#include <math.h> 
#include <string.h>
#include <ctime>

#define CombinLen 2      //SNP ��ϸ���
#define MaxSNPState 2    //SNP = 0,1 ... MaxSNPState
#define T_Value 0        //���ָ߸к͵͸е���ֵ
#define N_P     1000     //����Pֵ��ѭ������

const char path[100] = "D:\\Documents\\Visual Studio 2013\\Projects\\MDR\\data\\MDR-SampleData.txt";
const char k_path[100] = "D:\\Documents\\Visual Studio 2013\\Projects\\MDR\\data\\K.txt";

struct SNPState
{
	int SNP1;
	int SNP2;
	int State;              //0:��ȱ 1:�߸� 2:�͸�
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

	/*���Թ�ȷ���ļ��������������жϣ������ļ��Ĵ�С��ֵ*/
	char data_Name[data_Col][10];
	int data[data_Line][data_Col];

	// 1. ��ȡ����
	FILE *fp;
	fp = fopen(path, "r");

	for (int i = 0; i < data_Col; i++)    //��ȡ������
		fscanf(fp, "%s", data_Name[i]);

	for (int i = 0; i < data_Line; i++)   //��ȡ����
	{
		for (int j = 0; j < data_Col; j++)
		{
			fscanf(fp, "%d", &data[i][j]);
		}
	}
	fclose(fp);

	// 2. ѵ������
	//�������K���Ӽ�
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
	int K[data_Line][data_Col];                       //K�ǻ����Ӽ�������ݣ�˳��洢�˸����Ӽ� 
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

	//һ���Ӽ������Լ�������K-1���Ӽ���ѵ��������ѵ��������ѡѵ��׼ȷ����ߵ�N���Ա������
	//void KFoldCrossValid(int **data, const int row, const int col, SNPCombin *SNPData);
	//���룺ѵ�������ݣ�������������������
	//������Ա������  
	const int K_line = data_Line - data_Line / K_num;
	SNPCombin *SNPData_tmp, SNPData[K_num];
	SNPData_tmp = (struct SNPCombin*)malloc(sizeof(struct SNPCombin));              //�ṹ��ָ���ʼ��
	for (int k = 0; k < K_num; k++)
	{
		//����׼��
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

	//���㽻��һ����CVC
	//�����Ա�����ϸ���num
	//int num;
	//long int tmp1, tmp2 = 1;
	//tmp1 = 1;
	//for (int i = 0; i < CombinLen; i++)
	//{
	//tmp1 *=  data_Col - 1 - i;
	//tmp2 *=  CombinLen - i;
	//}
	//num = tmp1 / CombinLen; 

	//������Ա�����ϱ�ѡ��Ĵ���
	int *SNP_Count;
	SNP_Count = (int *)malloc((data_Col - 1) * (data_Col - 1) * sizeof(int));
	memset(SNP_Count, 0, (data_Col - 1) * (data_Col - 1) * sizeof(int));
	for (int i = 0; i < K_num; i++)
	{
		int tmp;
		tmp = SNPData[i].index_i * (data_Col - 1) + SNPData[i].index_j;
		SNP_Count[tmp]++;
	}
	//ѡ�����������,��CVC�����Ա������
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
	if (SNP_Count[max_SNP] < CVC_Threshold)             //ѡ�����ε�SNP��ϵĴ���Ӧ�ô���CVC��ֵ
	{
		printf("û�з����������Ա�����ϡ�\n");
		return;
	}
	else
		printf("�Ա������Ϊ��SNP%d  SNP%d\n", max_SNP1 + 1, max_SNP2 + 1);

	free(SNP_Count);

	//�ò����Ӽ�������ѡ�����Ա�����ϵ�׼ȷ��
	//void ExamSNP(int **data, const int row, const int col, SNPCombin SNPData, float ExamAccuracy);
	//���룺�����Ӽ��������Ӽ�������������������SNP��ϵĽṹ��
	//���������׼ȷ��
	const int Exam_line = data_Line / K_num;
	int Exam_data[Exam_line][data_Col];
	float ExamAccuracy[K_num] = { 0 };              //�洢����Ĳ���׼ȷ�ʡ�����ѡ�����ϲ����㣬׼ȷ��Ϊ0
	int exam_count = 0;
	for (int k = 0; k < K_num; k++)
	{
		int tmp = 0;
		if (SNPData[k].index_i == max_SNP1 && SNPData[k].index_j == max_SNP2)
		{
			// ����׼��
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

	//����Pֵ
	int P_count = 0;
	float P;
	for (int n = 0; n < N_P; n++)
	{
		for (int i = 0; i < data_Line; i++)
		{
			int tmp = rand() % 2;            //�Է�ֵ����û�
			data[i][data_Col - 1] = tmp;
		}

		//��ָ�����Ա�����Ͻ���ѵ���Ͳ���
		//���룺�Ա�����ϵ����������ݣ�����������
		//���������׼ȷ��
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

	printf("�������׼ȷ��Ϊ��%4.2f%% ƽ��׼ȷ��Ϊ��%4.2f%%\n", max_Accuracy * 100, avg_Accuracy * 100);
	printf("\nIf-Then Rules:\n");
	for (int i = 0; i < (MaxSNPState + 1)*(MaxSNPState + 1); i++)
	{
		printf("IF SNP1 = %d AND SNP2 = %d THEN CLASSIFY AS %d\n",
			SNPData[max_index].State[i].SNP1,
			SNPData[max_index].State[i].SNP2,
			SNPData[max_index].State[i].State - 1);
	}

	printf("\nPֵ:   %f\n", P);
	printf("\nCPU�����ʱ: %dms\n", time);
	printf("\n");
	return;
}

void KFoldCrossValid(int **data, const int row, const int col, SNPCombin *SNPData)
{
	//for (int i = 0; i < col; i++)
	//{
	//	printf("%d\t", *((int*)data + i));          //�ֶ�������ά�������ʶ��: *((int*)data + i*col + j) == data[i][j]
	//}



	//����ѵ��׼ȷ��
	float TrainAccuracy[200][3] = { 0 };                 //ѵ��׼ȷ�ʣ���ʱ������������壬����Ӧ��Ϊnum = 190 ��λ�㣬
	SNPState State[200][(MaxSNPState + 1)*(MaxSNPState + 1)];  //�洢�������
	int tmp1 = 0;                                             //tmp1:�ڼ����Ա���������ϵ�index
	for (int i = 0; i < col - 2; i++)
	{
		for (int j = i + 1; j < col - 1; j++)
		{
			float Case1Sum[(MaxSNPState + 1)*(MaxSNPState + 1)] = { 0 };
			float Case0Sum[(MaxSNPState + 1)*(MaxSNPState + 1)] = { 0 };                      // =9 �洢������ֵ�͡����շ�ֵ��

			int tmp2 = 0;                                    //tmp2:һ���Ա�����ϵĵڼ�����ʽ��index
			//����ÿһ��SNP�Ա�����ϵ�ѵ��׼ȷ��
			for (int k = 0; k < row; k++)
			{
				//�������ϵĲ�����ֵ�ͼ����շ�ֵ��
				tmp2 = *((int*)data + k*col + i) * (MaxSNPState + 1) + *((int*)data + k*col + j);
				if (*((int*)data + k*col + col - 1) == 1)
				{
					//Case1Sum[tmp2] += *((int*)data + k*col + col);       //���һ������
					Case1Sum[tmp2] ++;
				}
				else
				{
					//Case0Sum[tmp2] += *((int*)data + k*col + col);       //���һ������
					Case0Sum[tmp2] ++;
				}
			}

			//���ָ߸к͵͸У�������ʱ�����ò�ֵ����0����С��0����   case0:���� case1:����
			float Accuracy = 0;
			float AFSum = 0;
			for (int m = 0; m < (MaxSNPState + 1)*(MaxSNPState + 1); m++)
			{
				int SNP1 = m / (MaxSNPState + 1);
				int SNP2 = m % (MaxSNPState + 1);
				if (Case0Sum[m] == 0 && Case1Sum[m] == 0)             //��ȱ
				{
					State[tmp1][m].SNP1 = SNP1;
					State[tmp1][m].SNP2 = SNP2;
					State[tmp1][m].State = 0;
				}
				else
				{
					if (Case0Sum[m] - Case1Sum[m] >= T_Value)         //�߸�
					{
						State[tmp1][m].SNP1 = SNP1;
						State[tmp1][m].SNP2 = SNP2;
						State[tmp1][m].State = 1;
						Accuracy += Case0Sum[m];
					}
					else                                              //�͸�
					{
						State[tmp1][m].SNP1 = SNP1;
						State[tmp1][m].SNP2 = SNP2;
						State[tmp1][m].State = 2;
						Accuracy += Case1Sum[m];
					}
				}
				AFSum += Case0Sum[m] + Case1Sum[m];
			}

			//���ݲ�����ֵ�ͼ������ֵ�ͣ�����ѵ��׼ȷ��
			TrainAccuracy[tmp1][0] = Accuracy / AFSum;
			TrainAccuracy[tmp1][1] = i;         //�洢λ�ã�����������Ա������  
			TrainAccuracy[tmp1][2] = j;
			tmp1++;
		}
	}

	//ѡ��׼ȷ�������Ա������
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
	//J[MaxIndexi * (MaxSNPState + 1) + MaxIndexj] ++;          //��¼�Ա�����ϱ�ѡ�еĴ���
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
				if (*((int*)data + i*col + col) == SNPData.State[j].State - 1)          //����(0)��Ϊ�߸�(1),����(1)��Ϊ�͸�(2)
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
	float Case0Sum[(MaxSNPState + 1)*(MaxSNPState + 1)] = { 0 };                      // =9 �洢������ֵ�͡����շ�ֵ��

	for (int k = 0; k < row - row / k_num; k++)
	{
		//�������ϵĲ�����ֵ�ͼ����շ�ֵ��
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
	//���ָ߸к͵͸У�������ʱ�����ò�ֵ����0����С��0����   case0:���� case1:����
	//float Accuracy = 0;
	//float AFSum = 0;
	SNPState State[(MaxSNPState + 1)*(MaxSNPState + 1)];  //�洢�������
	for (int m = 0; m < (MaxSNPState + 1)*(MaxSNPState + 1); m++)
	{
		int SNP1 = m / (MaxSNPState + 1);
		int SNP2 = m % (MaxSNPState + 1);
		if (Case0Sum[m] == 0 && Case1Sum[m] == 0)             //��ȱ
		{
			State[m].SNP1 = SNP1;
			State[m].SNP2 = SNP2;
			State[m].State = 0;
		}
		else
		{
			if (Case0Sum[m] - Case1Sum[m] >= T_Value)         //�߸�
			{
				State[m].SNP1 = SNP1;
				State[m].SNP2 = SNP2;
				State[m].State = 1;
			}
			else                                              //�͸�
			{
				State[m].SNP1 = SNP1;
				State[m].SNP2 = SNP2;
				State[m].State = 2;
			}
		}
	}

	//����
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
				if (*((int*)data + i*col + col) == State[j].State - 1)          //����(0)��Ϊ�߸�(1),����(1)��Ϊ�͸�(2)
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