#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>

int num_of_genes;
int num_of_TFs;
int max_level;

int num_of_point_LD;
//int num_of_point_TD;

typedef struct Time_Course{
        char gene_ID[20];
        double *LD_exp;
        //double *TD_exp;
        int *level;
        //int TimePoint;
        //int network[10];
}Time_Course;

Time_Course *gene_exp_table;
Time_Course *TF_exp_table;

typedef struct Relation_table{
    char query_gene_ID[20];
    char target_gene_ID[20];
    int checked;
} R_table;

Relation_table *Pos_Coexp;

int num_of_pos_edge = 0;

float pos_cutoff_LD;
float neg_cutoff_LD;
float pos_no_cutoff = 0.5;
float neg_no_cutoff = -0.5;

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

    	//for(int k=0; k<10; k++) {
    	//	TF_exp_table[i].network[k] = 0;
    	//}
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

    	//for(int k=0; k<10; k++) {
    	//	gene_exp_table[i].network[k] = 0;
    	//}
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

void Read_TF_level (char *input) {
	FILE *fptr = fopen(input, "r");
	char GID[20];
	//char header1[20];
	//char header2[20];
	int LV;
	int num_of_lines = 0;

	max_level = 0;
	while(fscanf(fptr, "%s,", GID) != EOF) {
		fscanf(fptr,"\t%d", &LV);
		if (LV > max_level) max_level = LV;
		num_of_lines++;
	}
	
	rewind(fptr);

	for(int i=0; i<num_of_genes; i++) {
		gene_exp_table[i].level = new int[max_level+1];
		for(int j=0; j<max_level+1; j++) {
			gene_exp_table[i].level[j] = 0;
		}
	}
		
	for(int i=0; i<num_of_TFs; i++) {
		TF_exp_table[i].level = new int[1];
		TF_exp_table[i].level[0] = 0;
	}

	while(fscanf(fptr, "%s,", GID) != EOF) {
		fscanf(fptr,"\t%d", &LV);
		for(int i=0; i<num_of_TFs; i++) {
			if(strcmp(GID, TF_exp_table[i].gene_ID) == 0) {
				TF_exp_table[i].level[0] = LV;
			}
		}
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
 
void node_pair_generator_LD_or_TD() {
    
    FILE *fout1;
    FILE *fout2;
    
    double R_LD;

    //for(int i=0; i<num_of_TFs; i++) {
    //    for(int j=0; j<num_of_genes; j++) {
	//		for(int k=0; k<10; k++) {
	//			TF_exp_table[i].network[k] = 0;
				//gene_exp_table[j].network[k] = 0;
	//		}
	//	}
	//}

    for(int i=0; i<num_of_TFs; i++) {
        for(int j=0; j<num_of_genes; j++) {
            if(strcmp(TF_exp_table[i].gene_ID, gene_exp_table[j].gene_ID) != 0) {
                R_LD = r_calculator(i,j);
            
                if(R_LD >= pos_cutoff_LD) {
                    num_of_pos_edge++;
                }
            }
        }
    }

    Pos_Coexp = new R_table[num_of_pos_edge];
    
    for (int i=0; i<num_of_pos_edge; i++) {
        Pos_Coexp[i].checked = 0;
    }
    
    fout1 = fopen("Gene_list_in_each_level.csv","w");
    fout2 = fopen("Gene_level_matrix.csv","w");

    //fprintf(fout1, "TF gene ID, coexpression type, gene ID, PCC under C1\n");
    
    int index_pos = 0;
    
    
    
    for(int i=0; i<num_of_TFs; i++) {
        for(int j=0; j<num_of_genes; j++) {
            if(strcmp(TF_exp_table[i].gene_ID, gene_exp_table[j].gene_ID) != 0) {
                R_LD = r_calculator(i,j);
                    
                if(R_LD >= pos_cutoff_LD) {
                    //fprintf(fout1, "%s,%s,%lf\n", TF_exp_table[i].gene_ID, gene_exp_table[j].gene_ID, R_LD);
                    gene_exp_table[j].level[TF_exp_table[i].level[0]] = 1;
                    //strcpy(Pos_Coexp[index_pos].query_gene_ID, TF_exp_table[i].gene_ID);
                    //strcpy(Pos_Coexp[index_pos].target_gene_ID, gene_exp_table[j].gene_ID);
                    index_pos++;
                }
            }
        }
    }
    
    fprintf(fout1, "Gene ID, is TF, TO-GCN level with R = %0.2f\n", pos_cutoff_LD);
    
    for(int i=0; i<max_level; i++) {
    	for(int j=0; j<num_of_TFs; j++) {
    		if(TF_exp_table[j].level[0] == i+1) {
    			fprintf(fout1, "%s,%d,%d\n", TF_exp_table[j].gene_ID, 1, i+1);
    		}
    	}
    	for(int j=0; j<num_of_genes; j++) {
    		if(gene_exp_table[j].level[i+1] == 1) {
    			fprintf(fout1, "%s,%d,%d\n", gene_exp_table[j].gene_ID, 0, i+1);
    		}
    	}
    	//fprintf(fout1,"\n");
    }
    
    fprintf(fout2, "Gene ID");
    for(int i=1; i<max_level+1; i++) {
    	fprintf(fout2, ",level %d", i);
    }
    fprintf(fout2, "\n");
    
    for(int i=0; i<num_of_TFs; i++) {
    	fprintf(fout2, "%s", TF_exp_table[i].gene_ID);
    	for(int j=1; j<max_level+1; j++) {
    		if(TF_exp_table[i].level[0] == j) {
    			fprintf(fout2, ",1");
    		} else {
    			fprintf(fout2, ",0");
    		}
    	}
    	fprintf(fout2, "\n");
    }
    
    for(int i=0; i<num_of_genes; i++) {
    	fprintf(fout2, "%s", gene_exp_table[i].gene_ID);
    	for(int j=1; j<max_level+1; j++) {
    		fprintf(fout2, ",%d", gene_exp_table[i].level[j]);
    	}
    	fprintf(fout2, "\n");
    }
    
    fclose(fout1);
    fclose(fout2);
}

int main(int argc, char* argv[]) {
    char input_file1[100];   //TF gene list
    char input_file2[100];   //All gene list
    char input_file3[100];   //TF level
    
    if (argc != 6) {
        printf("\nUsage: GeneLevel No_of_time_points TF_gene_matrix non_TF_gene_matrix TF_level_table Cutoff\n\n");
        
    } else {
        
        num_of_point_LD = atoi(argv[1]);
        
        strcpy(input_file1, argv[2]);
        strcpy(input_file2, argv[3]);
        strcpy(input_file3, argv[4]);        
        pos_cutoff_LD = atof(argv[5]);

        FILE *fptr1 = fopen(input_file1, "r");
        FILE *fptr2 = fopen(input_file2, "r");
        FILE *fptr3 = fopen(input_file3, "r");
        
        if(fptr1 == NULL || fptr2 == NULL || fptr3 == NULL) {
            
            printf("\nCan't find input file. Please check the inupt file again!\n\n");
            
        } else {
        
     		fclose (fptr1);
    		fclose (fptr2);
    		fclose (fptr3);
        
            Read_Time_Course_Data_TFs(input_file1);
            Read_Time_Course_Data_genes(input_file2);
            Read_TF_level(input_file3);

            printf("No. of TFs: %d\n", num_of_TFs);
            printf("No. of non-TF genes: %d\n", num_of_genes);
            printf("No. of time points: %d\n", num_of_point_LD);
            //printf("No. of samples under Cond.2: %d\n", num_of_point_TD);
            printf("Cutoff: %1.2lf\n\n", pos_cutoff_LD);
            //printf("Generating tables of eight types of coexpression......\n");
        
            node_pair_generator_LD_or_TD();
        
            printf("Done!\n");
        }
    }
    
    
    return 0;
}
