#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>

int num_of_genes;
int num_of_TFs;

long long histogram_LD[201];

int num_of_point_LD;
//int num_of_point_TD;

char input_file_TF[100];   //TF gene list
char input_file_gene[100];   //All gene list

typedef struct Time_Course{
        char gene_ID[20];
        double *LD_exp;

}Time_Course;

Time_Course *gene_exp_table;
Time_Course *TF_exp_table;

void Read_Time_Course_Data_TFs (char *input) {
	FILE *fptr = fopen(input, "r");
    
	char GID[20];
	double LDE[num_of_point_LD];


	num_of_TFs = 0;
    while(fscanf(fptr,"%s", GID) != EOF) {
		for(int i=0; i<num_of_point_LD; i++) {
			fscanf(fptr,"\t%lf", &LDE[i]);
		}

		num_of_TFs++;
	}		

	rewind(fptr);
	//table initialization
	
	TF_exp_table = new Time_Course[num_of_TFs];
    
    for(int i=0; i<num_of_TFs; i++) {
        TF_exp_table[i].LD_exp = new double[num_of_point_LD];

    }
    
	for(int i=0; i<num_of_TFs; i++) {


        for(int j=0; j<num_of_point_LD; j++) {
            TF_exp_table[i].LD_exp[j] = 0;
        }

    }
    

	int index = 0;
    while(fscanf(fptr,"%s", GID) != EOF) {
		strcpy(TF_exp_table[index].gene_ID, GID);
		for(int i=0; i<num_of_point_LD; i++) {
			fscanf(fptr,"\t%lf", &LDE[i]);
			TF_exp_table[index].LD_exp[i] = LDE[i];
		}

		index++;
	}		
	
    fclose(fptr);    
}

void Read_Time_Course_Data_genes (char *input) {
	FILE *fptr = fopen(input, "r");

	char GID[20];
	double LDE[num_of_point_LD];

	
	num_of_genes = 0;
	while(fscanf(fptr,"%s", GID) != EOF) {
		for(int i=0; i<num_of_point_LD; i++) {
			fscanf(fptr,"\t%lf", &LDE[i]);
		}

		num_of_genes++;
	}	

	rewind(fptr);
	//table initialization
	gene_exp_table = new Time_Course[num_of_genes];
    
    for(int i=0; i<num_of_genes; i++) {
        gene_exp_table[i].LD_exp = new double[num_of_point_LD];

    }
    
	for(int i=0; i<num_of_genes; i++) {


        for(int j=0; j<num_of_point_LD; j++) {
            gene_exp_table[i].LD_exp[j] = 0;
        }

    }	

	int index = 0;
	while(fscanf(fptr,"%s", GID) != EOF) {    
		strcpy(gene_exp_table[index].gene_ID, GID);
		for(int i=0; i<num_of_point_LD; i++) {
			fscanf(fptr,"\t%lf", &LDE[i]);
			gene_exp_table[index].LD_exp[i] = LDE[i];
		}

		index++;
	}		
	
    fclose(fptr);
}

double r_calculator (int x, int y) {
    double N_LD = num_of_point_LD;
    double R, SUM_XY, SUM_X, SUM_Y, SUM_X2, SUM_Y2;
    double temp_XY, temp_X, temp_Y, temp_X2, temp_Y2;
    
    temp_XY = temp_X = temp_Y = temp_X2 = temp_Y2 = 0;

    for(int i=0; i<N_LD; i++) {
        temp_XY = temp_XY + (TF_exp_table[x].LD_exp[i] * gene_exp_table[y].LD_exp[i]);
        temp_X = temp_X + TF_exp_table[x].LD_exp[i];
        temp_Y = temp_Y + gene_exp_table[y].LD_exp[i];
        temp_X2 = temp_X2 + (TF_exp_table[x].LD_exp[i] * TF_exp_table[x].LD_exp[i]);
        temp_Y2 = temp_Y2 + (gene_exp_table[y].LD_exp[i] * gene_exp_table[y].LD_exp[i]);
    }
    
    SUM_XY = temp_XY;
    SUM_X = temp_X;
    SUM_Y = temp_Y;
    SUM_X2 = temp_X2;
    SUM_Y2 = temp_Y2;
    
    if (SUM_X == 0 || SUM_Y == 0) {
        R = 0;
    } else {
        R = (SUM_XY-(SUM_X)*(SUM_Y)/N_LD)/(sqrt((SUM_X2-(SUM_X)*(SUM_X)/N_LD)*(SUM_Y2-(SUM_Y)*(SUM_Y)/N_LD)));
    }
    return R;
}

void histogram_calculation (double x) {
    int y;
    y = (int)((x + 1) * 100);
    
    histogram_LD[y]++;

}

void node_pair_generator_LD_or_TD() {

    double R_LD;


    for(int i=0; i<num_of_TFs; i++) {
        for(int j=0; j<num_of_genes; j++) {
            if(strcmp(TF_exp_table[i].gene_ID, gene_exp_table[j].gene_ID) != 0) {
                R_LD = r_calculator(i,j);
                histogram_calculation(R_LD);

            }
        }
    }
}

void function_one () {
    FILE *fout1;
    double bin[201];
    double PDF_LD[201];
    double CDF_asc_LD[201];
    double CDF_desc_LD[201];


    long long sum_LD = 0;
    
    bin[0] = -1;
    
    for(int i=1; i<=200; i++) {
        bin[i] = bin[i-1] + 0.01;
    }
    
    for(int i=0; i<=200; i++) {
        sum_LD = sum_LD + histogram_LD[i];
    }
    
    
    for(int i=0; i<=200; i++) {
        PDF_LD[i] = (double) histogram_LD[i] / (double) sum_LD;
    }
    
    CDF_asc_LD[0] = PDF_LD[0];
    
    for(int i=1; i<=200; i++) {
        CDF_asc_LD[i] = CDF_asc_LD[i-1] + PDF_LD[i];
    }
    
    CDF_desc_LD[200] = PDF_LD[200];
    
    for(int i=199; i>=1; i--) {
        CDF_desc_LD[i] = CDF_desc_LD[i+1] + PDF_LD[i];
    }
    
    fout1 = fopen("PCC_histogram.tsv","w");

    fprintf(fout1, "BIN\t# in LD\tPDF\tCDF_asc\tCDF_desc\n");
    for(int i=0; i<=200; i++) {
        fprintf(fout1, "%lf\t%lld\t%lf\t%lf\t%lf\n", bin[i], histogram_LD[i], PDF_LD[i], CDF_asc_LD[i], CDF_desc_LD[i]);
    }
    
    
    fclose(fout1);
    
    int pos_cutoff_LD;
    int neg_cutoff_LD;
    int found_LD = 0;
    
    for (int i=200; i>=1; i--) {
        if(found_LD == 0) {
            if(CDF_desc_LD[i] > 0.05) {
                pos_cutoff_LD = i;
                found_LD = 1;
            }
        }

    }
    
    found_LD = 0;

    for (int i=0; i<=200; i++) {
        if(found_LD == 0) {
            if(CDF_asc_LD[i] > 0.05) {
                neg_cutoff_LD = i;
                found_LD = 1;
            }
        }
        
    }
    
    pos_cutoff_LD = pos_cutoff_LD + 2;
    
    printf("Recommended cutoff = %1.2lf\n", bin[pos_cutoff_LD]);
    //printf("Cutoff of negative coexpression = %1.2lf\n", bin[neg_cutoff_LD]);
}

 
int main(int argc, char* argv[]) {
    char input_file1[100];   //TF gene list
    char input_file2[100];   //All gene list
    
    if (argc != 3) {
        printf("\nUsage: Cutoff No_of_time_points TF_gene_matrix\n\n");
    } else {
        
        for (int i=0; i<=200; i++) {
            histogram_LD[i] = 0;
        }
        
        num_of_point_LD = atoi(argv[1]);
        
        strcpy(input_file1, argv[2]);
        strcpy(input_file2, argv[2]);


        FILE *fptr1 = fopen(input_file1, "r");
        
        if(fptr1 == NULL) {
        
            printf("\nCan't find the input file. Please check the inupt file again!\n\n");
            
            fclose (fptr1);
        
        } else {
        
            fclose (fptr1);
            
            Read_Time_Course_Data_TFs(input_file1);
            Read_Time_Course_Data_genes(input_file2);

            printf("No. of TFs: %d\n", num_of_genes);
            printf("No. of time points: %d\n\n", num_of_point_LD);
        
            //printf("Waiting for histogram generation and cutoff values calculation......\n\n");
            node_pair_generator_LD_or_TD();
            function_one();
        }
    }
    
    return 0;	
}
