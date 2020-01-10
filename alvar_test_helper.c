// Note: Preallocated short memory with memory optimisation
#include<string.h>
#include<stdlib.h>
#include<time.h>
#include<math.h>
#include<omp.h>
#include<stdio.h>

float* FOAV(float* x, long len, int seq_len, int threadNum)
{
    float *av = malloc((len/2+1)*sizeof(float));
	int *p;
	switch(threadNum)
	{
		case 1: 
		{
			p = (int [2]){0,seq_len/2-1};	
			break;
		}
		case 2: 
		{
			p = (int [3]){0,seq_len*2.3/16,seq_len/2-1};
			break;
		}
		case 4:
		{
			p = (int [5]){0,seq_len*1.05/16,seq_len*2.3/16,seq_len*4/16,seq_len/2-1};
			break;
		}
		default: 
			return NULL;
	}
	float *short1[threadNum];
	float *short2[threadNum];
	for(int k = 0; k < threadNum; k++)
	{
		short1[k] = malloc(sizeof(float)*(len-2*p[k]+1));
		short2[k] = malloc(sizeof(float)*(len-2*p[k]+1));
		memset(short1[k], 0, (len-2*p[k]+1));
		memset(short2[k], 0, (len-2*p[k]+1));
	}
	clock_t end = clock();
	av[0]=0.0f;	
	#pragma omp parallel for
	for(int k = 0 ; k < threadNum; k++)
	{
        printf("AV started in %d\n",k);
		int limit = len-2*p[k], count = 4+2*p[k];
	    for(int j = 1; j <= p[k]+1;j++)// Initilise first mean
		{
			short1[k][1] += x[j];
	        short2[k][1] += x[p[k]+j+1];
        }	
        for(int j = 2; j <= len-2*p[k];j++)// Initilise successive means
		{
			short1[k][j] = short1[k][j-1] - x[j-1] + x[p[k]+j];
	       	short2[k][j] = short2[k][j-1] - x[p[k]+j] + x[2*p[k]+j+1];
		}
	    for(int n = p[k]+1 ; n <= p[k+1]; n ++, count += 2, limit -= 2)// calculate FOAV
		{	
			float denom = (float) (2*(len-2*n+1))*pow(n,2), sum = 0.0f;
      		for(int j = 1; j < limit; j++)
			{
					sum += (short1[k][j]-short2[k][j]) * (short1[k][j]-short2[k][j]);
                	short1[k][j] = short1[k][j] + x[n+j];
		        	short2[k][j] = short2[k][j+1] + x[count+j-1];
	        }
	        av[n] = sum/denom;
        }
    	clock_t end1 = clock();
		double time_spent4 = (double)(end1 - end) / CLOCKS_PER_SEC;
  		printf("seg %d takes %f ! \n",k,time_spent4);
        printf(" 1 completed by thread %d \n",omp_get_thread_num());// Time each thread
 	}
    for(int j = 0; j<threadNum; j++)
	{
		free(short1[j]);
		free(short2[j]);
	}
     return av;
}
int main()
{
	FILE *ptr_file, *fp;
	int j;
    const int len=3000;
    const int threadNum=1;
	char *bufx = (char *)malloc(10*sizeof(char));
    float *x;
	x = (float *)malloc((len+1)*sizeof(float));
	ptr_file = fopen("large.txt", "r");
	fp = fopen("allvar_large.txt", "w");
	if (!ptr_file)
		return 0;
	j=0;
	clock_t begin = clock();
	fscanf(ptr_file, "%s",bufx);
	x[j++]=0.0;
	x[j++]=atof(bufx);
	
	while (j<(len+1))
	{
		fscanf(ptr_file,"%s",bufx);
		x[j]=atof(bufx);
		j++;
	}
	printf("function called \n");
    float *av=FOAV( x, len, 20000, threadNum);
 	
	printf("\nAV ends at len = %d to %d \n",len,len/2-1);
	for(j = 0; j<10; j++)
	{
		printf("%f\n",av[j]);
	}    
    j = 0;
	while(j < len/2)
	{
		fprintf(fp, "%f \t",av[j++]);
	}
	clock_t end_total = clock();
	float time_spent2 = (float)(end_total - begin) / CLOCKS_PER_SEC;
    printf(" Total time  %f ! \n",time_spent2);
	fclose(fp);
	fclose(ptr_file);
	free(bufx);
	free(x);
	free(av);
    printf("files closed");
    return 0;
}
