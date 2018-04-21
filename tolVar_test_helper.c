// to calculate Total Variance for a given gyroscope data
// extracting from 54col_set1.txt
/* Note: len has to adjusted for each of the text files
*/
#include<string.h>
#include<stdlib.h>
#include<time.h>
#include<math.h>
#include<omp.h>
#include<stdio.h>
float* tvcalc(float *x_new, int len, int seq_len, int threadNum)
{
	
	const int N 	=	len*3-2;
	float *x		=	(float *)malloc((N+1)*sizeof(float));
	float *tv		=	(float *)malloc((len+1)*sizeof(float));
    
	for(int k=1;k<len;k++)
	{
		x[k]		=	x_new[len-k];		
		x[len-1+k]		=	x_new[k];
		x[2*len-1+k]	 =	x_new[len+1-k]; 
	}	
	x[2*len-1] = x_new[len]; 
	printf("TV started 3 \n");
	//int *p;
	int *p;
	switch(threadNum)
	{
		case 1: 
		{
			p = (int [2]){0, seq_len/2-1};	
			break;
		}
		case 2: 
		{
			p = (int [3]){0, seq_len*2.3/16, seq_len/2-1};
			break;
		}
		case 4:
		{
			p = (int [5]){0, seq_len*1.05/16, seq_len*2.3/16, seq_len*4/16, seq_len/2-1};
			break;
		}
		default: 
			return NULL;
	}
	for(int j=0;j<=threadNum;j++)
		printf("%d\t",p[j]);

//	printf("\n");

	float *short1[threadNum];
	float *short2[threadNum];
	for(int k = 0; k < threadNum; k++)
	{
		
		short1[k] = malloc(sizeof(float)*(N-2*p[k]+1));
		short2[k] = malloc(sizeof(float)*(N-2*p[k]+1));
		memset(short1[k], 0, (N-2*p[k]+1));
		memset(short2[k], 0, (N-2*p[k]+1));
	}
	printf("TV started 5 \n");
	clock_t end = clock();
	tv[0]=0.0f;	
	#pragma omp parallel for
	for(int k = 0 ; k < threadNum; k++)
	{
        printf("TV started in %d\n",k);
		int limit = N-2*p[k], count = 4+2*p[k];
	    for(int j = 1; j <= p[k]+1;j++)// Initilise first mean
		{
			short1[k][1] += x[j];
	        short2[k][1] += x[p[k]+j+1];
        }	
        for(int j = 2; j <= N-2*p[k];j++)// Initilise successive means
		{
			short1[k][j] = short1[k][j-1] - x[j-1] + x[p[k]+j];
	       	short2[k][j] = short2[k][j-1] - x[p[k]+j] + x[2*p[k]+j+1];
		}
	    for(int n = p[k]+1 ; n <= p[k+1]; n ++, count += 2, limit -= 2)// claculate FOAV
		{	
			float denom	=	(float) (2*(3*len-2*n-1))*pow(n,2);
			float sum = 0.0f;
      		for(int j = 1; j < limit; j++)
			{
					sum += (short1[k][j]-short2[k][j]) * (short1[k][j]-short2[k][j]);
                	short1[k][j] = short1[k][j] + x[n+j];
		        	short2[k][j] = short2[k][j+1] + x[count+j-1];
	        }
	        tv[n] = sum/denom;
        }
    	clock_t end1 = clock();
		double time_spent4 = (double)(end1 - end) / CLOCKS_PER_SEC;
  		printf("seg %d takes %f ! \n",k,time_spent4);
        printf(" 1 completed by thread %d \n",omp_get_thread_num());// Time each thread
 	}
	printf("TV started 5 \n");
    for(int j = 0; j<threadNum; j++)
	{
		free(short1[j]);
		printf("TV 1 \n");
		free(short2[j]);
		printf("TV 2 %d \n",j);
	}
	printf("TV started 6 \n");
	return tv;
}
int main()
{
	FILE *ptr_file,*fp;
	const int len=30000;//238256;
	char bufx[10];
	ptr_file = fopen("large.txt","r");
	if (!ptr_file) return 0;

	float *x_new;
	x_new = (float *)malloc((len+1)*sizeof(float));

	fscanf(ptr_file, "%s",bufx);
	int j=1;	x_new[j++] = 0;
	while (j<(len+1)){
		fscanf(ptr_file, "%s",bufx);
		x_new[j++] = atof(bufx);
	}
	fclose(ptr_file);

// Time segment 1 
	clock_t start = clock();
/**/  

	float *tv = tvcalc(x_new, len, 20000,2);

// Time segment 2
	clock_t end = clock();
        double time_spent2 = (double)(end - start) / CLOCKS_PER_SEC;
	printf(" Total variance done in %f ! \n",time_spent2);
/**/

	printf("TV ends \n");
	for(int j=0;j<10;j++){
		printf	("%f\n",tv[j]);
	}    
	fp = fopen("TolVar_op.txt", "w+");
	j=0;
	printf("%d",len);
	while(j<len){	
		fprintf(fp, "%f \n",tv[j++]);
	}
	printf("\n fp");
	fclose(fp);
	printf("\n xnew \n");
	free(x_new);
	free(tv);
	printf("tv");
	
	printf("files closed");
    return 0;
}