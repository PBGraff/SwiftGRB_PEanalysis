// gcc-mp-4.8 -O3 -o test kstest.c test.c 

#include <string.h>
#include "kstest.h"

int main(int argc, char **argv)
{
	//SkyNetKStest(argv[1]);

	char filename[200];
	int i;
	/*
	char nhid[16][20] = {"25","25","25-25","25-25","50","50","50-50","50-50","100","100","100-30","100-30","100-50","100-50","100-100","100-100"};
	char act[16][5] = {"10","30","110","330","10","30","110","330","10","30","110","330","110","330","110","330"};
	for ( i=0; i<16; i++ )
	{
		strcpy(filename,"../../SwiftGRBpipe/NNpreds/priorsample2_CVall_nhid-");
		strcat(filename,nhid[i]);
		strcat(filename,"_act");
		strcat(filename,act[i]);
		strcat(filename,"_blind_pred.txt");
		printf("%s\t%s\t",nhid[i],act[i]);
		SkyNetKStest(filename);
	}
	*/
	printf("\n");
	char SKLmodels[4][30] = {"RandomForest","AdaBoost","SupportVectorMachine","NeuralNetwork"};
	for ( i=0; i<4; i++ )
	{
		strcpy(filename,"../models/");
		strcat(filename,SKLmodels[i]);
		strcat(filename,"_predictionsKS.txt");
		printf("%s\t",SKLmodels[i]);
		if (i==2) printf("\t");
		SciKitLearnKStest(filename);
	}

	return 0;
}

