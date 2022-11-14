#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include "sndfile.h"
#include <math.h>

#define N_interpol 7
#define N_decimacao 7
#define N_aliasing 7
#define N_filtragens_interpol 3
#define N_filtragens_decimacao 3
#define N_distorcao_max 10000

void Interpol(double* amostra_interpol_in, double* amostra_interpol_out, double *entrada_interpol,double *saida_interpol,double *a_interpol, double *b_interpol)
{    

    double soma_entrada;
    double soma_saida;
    double entrada_multiplicada[N_interpol+1];
    double saida_multiplicada[N_interpol];

    amostra_interpol_out[0]=amostra_interpol_in[0];
    
    for(int k=0; k<N_filtragens_interpol; k++)
    {
        soma_entrada=0;
        soma_saida=0;
        //Reordenando Vetor de entrada para próxima saída
        for(int i=N_interpol; i>0; i--)
        {
            entrada_interpol[i+k*(N_interpol+1)]=entrada_interpol[i-1+k*(N_interpol+1)];
        }
        entrada_interpol[k*(N_interpol+1)]=amostra_interpol_out[0];

        //Calculo da saida
        for(int i=0; i<N_interpol+1; i++)
        {
            entrada_multiplicada[i] = entrada_interpol[i+k*(N_interpol+1)]*b_interpol[i];
            soma_entrada=soma_entrada+entrada_multiplicada[i];
        }

        for(int i=0; i<N_interpol; i++)
        {
            saida_multiplicada[i] = saida_interpol[i+k*(N_interpol)]*a_interpol[i];
            soma_saida=soma_saida+saida_multiplicada[i];
        }

        amostra_interpol_out[0]=soma_entrada-soma_saida;

        //Reordenando Vetor de saida para próxima saída
        for(int i=N_interpol-1; i>0; i--)
        {
            saida_interpol[i+k*(N_interpol)]=saida_interpol[i-1+k*(N_interpol)];
        }
        
        saida_interpol[k*(N_interpol)]=amostra_interpol_out[0];
    }

    //Controle dos limites da saída entre +1 e -1
    if(amostra_interpol_out[0]>1)
    {
        amostra_interpol_out[0]=1;
    }
    else if(amostra_interpol_out[0]<-1)
    {
        amostra_interpol_out[0]=-1;
    }

}

void Decimacao(double* amostra_decimacao_in, double* amostra_decimacao_out, double *entrada_decimacao,double *saida_decimacao,double *a_decimacao, double *b_decimacao)
{
    double soma_entrada;
    double soma_saida;
    double entrada_multiplicada[N_decimacao+1];
    double saida_multiplicada[N_decimacao];

    amostra_decimacao_out[0] = amostra_decimacao_in[0];

    for(int k=0; k<N_filtragens_decimacao; k++)
    {
        soma_entrada=0;
        soma_saida=0;

        //Reordenando Vetor de entrada para próxima saída
        for(int i=N_decimacao; i>0; i--)
        {
            entrada_decimacao[i+k*(N_decimacao+1)]=entrada_decimacao[i-1+k*(N_decimacao+1)];
        }
        entrada_decimacao[k*(N_decimacao+1)]=amostra_decimacao_out[0];

        //Calculo da saida
        for(int i=0; i<N_decimacao+1; i++)
        {
            entrada_multiplicada[i] = entrada_decimacao[i+k*(N_decimacao+1)]*b_decimacao[i];
            soma_entrada=soma_entrada+entrada_multiplicada[i];
        }

        for(int i=0; i<N_decimacao; i++)
        {
            saida_multiplicada[i] = saida_decimacao[i+k*(N_decimacao)]*a_decimacao[i];
            soma_saida=soma_saida+saida_multiplicada[i];
        }

        
        amostra_decimacao_out[0] = soma_entrada-soma_saida;

        //Reordenando Vetor de saida para próxima saída
        for(int i=N_decimacao-1; i>0; i--)
        {
            saida_decimacao[i+k*(N_decimacao)]=saida_decimacao[i-1+k*(N_decimacao)];
        }

        saida_decimacao[k*(N_decimacao)]=amostra_decimacao_out[0];

    }
}

void Distorcao_filtrada(double* amostra_aliasing_in, double* amostra_aliasing_out, double *entrada_aliasing,double *saida_aliasing,double *a_aliasing, double *b_aliasing, double grau_distorcao)
{    

    double soma_entrada;
    double soma_saida;
    double entrada_multiplicada[N_aliasing+1];
    double saida_multiplicada[N_aliasing];

    //Reordenando Vetor de entrada para próxima saída
    soma_entrada=0;
    soma_saida=0;
    for(int i=N_aliasing; i>0; i--)
    {
        entrada_aliasing[i]=entrada_aliasing[i-1];
    }
    entrada_aliasing[0]=atan(grau_distorcao*amostra_aliasing_in[0])/atan(grau_distorcao);

    //printf("dps reoordenacao entrada\n");

    //Calculo da saida
    for(int i=0; i<N_aliasing+1; i++)
    {
        entrada_multiplicada[i] = entrada_aliasing[i]*b_aliasing[i];
        soma_entrada=soma_entrada+entrada_multiplicada[i];
    }

    //printf("dps calc entrada mult\n");

    for(int i=0; i<N_aliasing; i++)
    {
        saida_multiplicada[i] = saida_aliasing[i]*a_aliasing[i];
        soma_saida=soma_saida+saida_multiplicada[i];
    }

    //printf("dps calc saida mult\n");

    amostra_aliasing_out[0]=soma_entrada-soma_saida;

   
    //Reordenando Vetor de saida para próxima saída
    for(int i=N_aliasing-1; i>0; i--)
    {
        saida_aliasing[i]=saida_aliasing[i-1];
    }

    //printf("dps reoordenacao saida\n");

    saida_aliasing[0]=amostra_aliasing_out[0];
}



int main(int argc, char *argv[])
{
    if (argc != 2)
    {
        printf("./executavel [N_repeticoes]");
        return 1;
    }
    
    //Definições do arquivo de leitura e saida
    SNDFILE *file_in,*file_out ;   //Arquivo de entrada e saída
    SF_INFO sfinfo_in,sfinfo_out ; //Arquivo de informações de entrada e saída
    
    
    //Definições do Projeto
    const int L = 30; //Fator de interpolação
    const double grau_distorcao=0.1; //Grau de distorcao

    //Variaveis da CPU
    double *amostra_in;
    double *aux;
    double *amostra_out;
    double *entrada_interpol;
    double *saida_interpol;
    double *entrada_decimacao;
    double *saida_decimacao;
    double *entrada_aliasing;
    double *saida_aliasing;

    //Alocação das Variaveis da CPU
    amostra_in = (double*)malloc(sizeof(double) * 1);
    aux = (double*)malloc(sizeof(double) * 1);
    amostra_out = (double*)malloc(sizeof(double) * 1);
    entrada_interpol = (double*)malloc(sizeof(double) * (N_filtragens_interpol)*(N_interpol+1));
    saida_interpol = (double*)malloc(sizeof(double) * (N_filtragens_interpol)*N_interpol);
    entrada_decimacao = (double*)malloc(sizeof(double) * (N_filtragens_decimacao)*(N_decimacao+1));
    saida_decimacao = (double*)malloc(sizeof(double) * (N_filtragens_decimacao)*N_decimacao);
    entrada_aliasing = (double*)malloc(sizeof(double) * (N_distorcao_max*(N_aliasing+1)));
    saida_aliasing = (double*)malloc(sizeof(double) * (N_distorcao_max*N_aliasing));


    memset(entrada_interpol, 0, (N_filtragens_interpol)*(N_interpol+1)*sizeof(double));
    memset(entrada_decimacao, 0, (N_filtragens_decimacao)*(N_decimacao+1)*sizeof(double));
    memset(entrada_aliasing, 0, N_distorcao_max*(N_aliasing+1)*sizeof(double));

    memset(saida_interpol, 0, (N_filtragens_interpol)*N_interpol*sizeof(double));
    memset(saida_decimacao, 0, (N_filtragens_decimacao)*N_decimacao*sizeof(double));
    memset(saida_aliasing, 0, N_distorcao_max*N_aliasing*sizeof(double));


    //Definição dos coeficientes do filtro na CPU
    double a_interpol[N_interpol] = {
        #include "a_interpol_30.txt"
    };
    double b_interpol[N_interpol+1] = {
        #include "b_interpol_30.txt"
    };
    double a_decimacao[N_decimacao] = {
        #include "a_decimacao_30.txt"
    };
    double b_decimacao[N_decimacao+1] = {
        #include "b_decimacao_30.txt"
    };
    double a_aliasing[N_aliasing] = {
        #include "a_aliasing_30.txt"
    };
    double b_aliasing[N_aliasing+1] = {
        #include "b_aliasing_30.txt"
    };
    
    clock_t inicio_total;
	clock_t fim_total;

    //Leitura de amostras e processamento
    double time = 0;
    double soma_total=0;
    int num_iteracoes=0;
    int read_count=1;

    FILE *fp = fopen("time_cpu_L=30_5k", "w");
    FILE *fp_amostras = fopen("time_cpu_L=30_5k_amostras", "w");
    fprintf(fp, "N_distorcao, Tempo");
    fprintf(fp_amostras, "N_distorcao, Tempo");
    
    int nro_repeticoes_de_amostra = atoi(argv[1]);

    for(int N_distorcao = 144; N_distorcao<=5040; N_distorcao=N_distorcao+144)
    {
        soma_total=0;

        for(int z=0; z<nro_repeticoes_de_amostra;z++)
        {
            num_iteracoes=0;

            memset(entrada_interpol, 0, (N_filtragens_interpol)*(N_interpol+1)*sizeof(double));
            memset(entrada_decimacao, 0, (N_filtragens_decimacao)*(N_decimacao+1)*sizeof(double));
            memset(entrada_aliasing, 0, N_distorcao_max*(N_aliasing+1)*sizeof(double));

            memset(saida_interpol, 0, (N_filtragens_interpol)*N_interpol*sizeof(double));
            memset(saida_decimacao, 0, (N_filtragens_decimacao)*N_decimacao*sizeof(double));
            memset(saida_aliasing, 0, N_distorcao_max*N_aliasing*sizeof(double));

            sfinfo_in.format = 0;          //Documentação do libsnd manda fazer isso para arquivos de leitura
            file_in = sf_open ("sweep_48khz_5s.wav", SFM_READ, &sfinfo_in); //Determina o arquivo de entrada
            sfinfo_out=sfinfo_in;
            file_out = sf_open ("sweep_48khz_time_cpu.wav", SFM_WRITE, &sfinfo_out); //Determina o arquivo de saída
            sf_command (file_out, SFC_SET_CLIPPING, NULL, SF_TRUE) ;

            inicio_total = clock();  

            while (num_iteracoes<5000)
            {   
                read_count = (int) sf_read_double (file_in, amostra_in, 1);
                
                //Interpolacao i
                Interpol(&amostra_in[0], &amostra_out[0], entrada_interpol, saida_interpol, a_interpol, b_interpol);


                for(int j=0; j<N_distorcao; j++)
                {
                    Distorcao_filtrada(&amostra_out[0],&amostra_out[0],&entrada_aliasing[(N_aliasing+1)*j],&saida_aliasing[(N_aliasing)*j],a_aliasing,b_aliasing,grau_distorcao);
                }

                //Decimacao i
                Decimacao(&amostra_out[0], &amostra_out[0], entrada_decimacao, saida_decimacao, a_decimacao, b_decimacao);

                aux[0] = amostra_out[0];

                for(int i=1; i<L; i++)
                {
                    //Interpolacao i
                    Interpol(&amostra_in[0], &amostra_out[0], entrada_interpol, saida_interpol, a_interpol, b_interpol);

                    for(int j=0; j<N_distorcao; j++)
                    {
                        Distorcao_filtrada(&amostra_out[0],&amostra_out[0],&entrada_aliasing[(N_aliasing+1)*j],&saida_aliasing[(N_aliasing)*j],a_aliasing,b_aliasing,grau_distorcao);
                    }

                    //Decimacao i
                    Decimacao(&amostra_out[0], &amostra_out[0], entrada_decimacao, saida_decimacao, a_decimacao, b_decimacao);
                }

                // Escrevendo no arquivo de saida
                sf_write_double (file_out, &aux[0], read_count) ;
                num_iteracoes++;
            }

            fim_total = clock();

            time = (((double)(fim_total - inicio_total)) / CLOCKS_PER_SEC);
            
            fprintf(fp, "%d, %f, %f \n",N_distorcao,time);
            printf("N_distorcoes: %d Tempo [s]: %f \n", N_distorcao, time);

            soma_total+=time;

            sf_close(file_in);
            sf_close(file_out);
        }

        soma_total=soma_total/nro_repeticoes_de_amostra;
        
        fprintf(fp, "%d, %f, %f \n",N_distorcao,soma_total);
    }
    


    //Liberando a memória alocada
    free(amostra_in);
    free(aux);
    free(amostra_out);
    fclose(fp);

    return 0;
}