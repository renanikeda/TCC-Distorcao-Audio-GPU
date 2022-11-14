#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include "sndfile.h"
#include <math.h>

#define N_interpol 6
#define N_decimacao 6
#define N_aliasing 7
#define N_estagios 2000

void Interpol(double* amostra_interpol_in, double* amostra_interpol_out, double *entrada_interpol,double *saida_interpol,double *a_interpol, double *b_interpol)
{    

    double soma_entrada;
    double soma_saida;
    double entrada_multiplicada[N_interpol+1];
    double saida_multiplicada[N_interpol];

    soma_entrada=0;
    soma_saida=0;
    //Reordenando Vetor de entrada para próxima saída
    for(int i=N_interpol; i>0; i--)
    {
        entrada_interpol[i]=entrada_interpol[i-1];
    }
    entrada_interpol[0]=amostra_interpol_in[0];

    //Calculo da saida
    for(int i=0; i<N_interpol+1; i++)
    {
        entrada_multiplicada[i] = entrada_interpol[i]*b_interpol[i];
        soma_entrada=soma_entrada+entrada_multiplicada[i];
    }

    for(int i=0; i<N_interpol; i++)
    {
        saida_multiplicada[i] = saida_interpol[i]*a_interpol[i];
        soma_saida=soma_saida+saida_multiplicada[i];
    }

    amostra_interpol_out[0]=soma_entrada-soma_saida;

    //Reordenando Vetor de saida para próxima saída
    for(int i=N_interpol-1; i>0; i--)
    {
        saida_interpol[i]=saida_interpol[i-1];
    }
    
    saida_interpol[0]=amostra_interpol_out[0];

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

    soma_entrada=0;
    soma_saida=0;

    //Reordenando Vetor de entrada para próxima saída
    for(int i=N_decimacao; i>0; i--)
    {
        entrada_decimacao[i]=entrada_decimacao[i-1];
    }
    entrada_decimacao[0]=amostra_decimacao_in[0];

    //Calculo da saida
    for(int i=0; i<N_decimacao+1; i++)
    {
        entrada_multiplicada[i] = entrada_decimacao[i]*b_decimacao[i];
        soma_entrada=soma_entrada+entrada_multiplicada[i];
    }

    for(int i=0; i<N_decimacao; i++)
    {
        saida_multiplicada[i] = saida_decimacao[i]*a_decimacao[i];
        soma_saida=soma_saida+saida_multiplicada[i];
    }

    
    amostra_decimacao_out[0] = soma_entrada-soma_saida;

    //Reordenando Vetor de saida para próxima saída
    for(int i=N_decimacao-1; i>0; i--)
    {
        saida_decimacao[i]=saida_decimacao[i-1];
    }

    saida_decimacao[0]=amostra_decimacao_out[0];
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



int main()
{
    //Definições do arquivo de leitura e saida
    SNDFILE *file_in,*file_out ;   //Arquivo de entrada e saída
    SF_INFO sfinfo_in,sfinfo_out ; //Arquivo de informações de entrada e saída
    sfinfo_in.format = 0;          //Documentação do libsnd manda fazer isso para arquivos de leitura
    file_in = sf_open ("sweep_48khz_5s.wav", SFM_READ, &sfinfo_in); //Determina o arquivo de entrada
    sfinfo_out=sfinfo_in;
    file_out = sf_open ("sweep_48khz_5s_2000_0.045.wav", SFM_WRITE, &sfinfo_out); //Determina o arquivo de saída
    sf_command (file_out, SFC_SET_CLIPPING, NULL, SF_TRUE) ;
    printf("Informações sobre o arquivo de entrada: \n");
    printf("Taxa de amostragem = %d , Frames = % d , Canais = % d \n" , (int) sfinfo_in.samplerate,(int) sfinfo_in.frames,(int)        	  sfinfo_in.channels);                                                              //Mostra algumas caratersiticas do arquivo de entrada 
    
    //Definições do Projeto
    const int L = 30; //Fator de interpolação
    const double grau_distorcao=0.045; //Grau de distorcao

    //Variaveis da CPU
    double *amostra_in;
    double *amostra_out;
    double *entrada_interpol;
    double *saida_interpol;
    double *entrada_decimacao;
    double *saida_decimacao;
    double *entrada_aliasing;
    double *saida_aliasing;

    //Alocação das Variaveis da CPU
    amostra_in = (double*)malloc(sizeof(double) * 1);
    amostra_out = (double*)malloc(sizeof(double) * 1);
    entrada_interpol = (double*)malloc(sizeof(double) * (N_interpol+1));
    saida_interpol = (double*)malloc(sizeof(double) * N_interpol);
    entrada_decimacao = (double*)malloc(sizeof(double) * (N_decimacao+1));
    saida_decimacao = (double*)malloc(sizeof(double) * N_decimacao);
    entrada_aliasing = (double*)malloc(sizeof(double) * (N_estagios*(N_aliasing+1)));
    saida_aliasing = (double*)malloc(sizeof(double) * (N_estagios*N_aliasing));


    memset(entrada_interpol, 0, (N_interpol+1)*sizeof(double));
    memset(entrada_decimacao, 0, (N_decimacao+1)*sizeof(double));
    memset(entrada_aliasing, 0, N_estagios*(N_aliasing+1)*sizeof(double));

    memset(saida_interpol, 0, N_interpol*sizeof(double));
    memset(saida_decimacao, 0, N_decimacao*sizeof(double));
    memset(saida_aliasing, 0, N_estagios*N_aliasing*sizeof(double));


    //Definição dos coeficientes do filtro na CPU
    double a_interpol_host[N_interpol] = {
        #include "a_interpol_30_final.txt"
    };
    double b_interpol_host[N_interpol_mais_1] = {
        #include "b_interpol_30_final.txt"
    };
    double a_decimacao_host[N_decimacao] = {
        #include "a_decimacao_30_final.txt"
    };
    double b_decimacao_host[N_decimacao_mais_1] = {
        #include "b_decimacao_30_final.txt"
    };
    double a_aliasing_host[N_aliasing] = {
        #include "a_aliasing_30_final.txt"
    };
    double b_aliasing_host[N_aliasing_mais_1] = {
        #include "b_aliasing_30_final.txt"
    };
    
    clock_t inicio_total;
	clock_t fim_total;

    //Leitura de amostras e processamento
    int read_count = 1;
    int num_iteracoes = 0;
    int pct=0;
    double soma_total = 0;

    printf("Começando o loop...\n");

    while (read_count)
    {
        inicio_total = clock();       
        read_count = (int) sf_read_double (file_in, amostra_in, 1);
        
        //Interpolacao i
        Interpol(&amostra_in[0], &amostra_out[0], entrada_interpol, saida_interpol, a_interpol, b_interpol);


        for(int j=0; j<N_estagios; j++)
        {
            Distorcao_filtrada(&amostra_out[0],&amostra_out[0],&entrada_aliasing[(N_aliasing+1)*j],&saida_aliasing[(N_aliasing)*j],a_aliasing,b_aliasing,grau_distorcao);
        }

        //Decimacao i
        Decimacao(&amostra_out[0], &amostra_out[0], entrada_decimacao, saida_decimacao, a_decimacao, b_decimacao);

        // Escrevendo no arquivo de saida
        sf_write_double (file_out, &amostra_out[0], read_count) ;

        for(int i=1; i<L; i++)
        {
            //Interpolacao i
            Interpol(&amostra_in[0], &amostra_out[0], entrada_interpol, saida_interpol, a_interpol, b_interpol);

            for(int j=0; j<N_estagios; j++)
            {
                Distorcao_filtrada(&amostra_out[0],&amostra_out[0],&entrada_aliasing[(N_aliasing+1)*j],&saida_aliasing[(N_aliasing)*j],a_aliasing,b_aliasing,grau_distorcao);
            }

            //Decimacao i
            Decimacao(&amostra_out[0], &amostra_out[0], entrada_decimacao, saida_decimacao, a_decimacao, b_decimacao);
        }

        fim_total = clock();

        soma_total += ((double)(fim_total - inicio_total)) / CLOCKS_PER_SEC;
        
        if(num_iteracoes % 17640 == 0)
        {
            printf("Porcentagem do processamento = %d%c \n", pct,37);
            pct+=10;
        }
        // if(num_iteracoes > 9)
        // {
        //    break;
        // }
        num_iteracoes++;
    }

    printf("Acabou\n");

    printf("Nro de amostras processados.....: %d\n", num_iteracoes);
    printf("\n");
    printf("Tempo Total \n"); 
    printf("Tempo total de processamento [s]...: %f \n", soma_total);
    // printf("Tempo total de interpolacao [ms]...: %f \n", soma_interpolacao);
    // printf("Tempo total de decimacao [ms]...: %f \n", soma_decimacao);
    // printf("Tempo total de transferencia de memoria H/D [ms]...: %f \n", soma_thd);
    // printf("Tempo total de transferencia de memoria D/H [ms]...: %f \n", soma_tdh);

    printf("\n");
    printf("Tempo Medio \n");
    printf("Media do tempo de processamento [s]...: %f \n", soma_total/num_iteracoes);
    // printf("Media do tempo total de interpolacao [ms]...: %f \n", soma_interpolacao/num_iteracoes);
    // printf("Media do tempo total de decimacao [ms]...: %f \n", soma_decimacao/num_iteracoes);
    // printf("Media do tempo total de transferencia de memoria H/D [ms]...: %f \n", L*soma_thd/num_iteracoes);
    // printf("Media do tempo total de transferencia de memoria D/H [ms]...: %f \n", L*soma_tdh/num_iteracoes);

    printf("\n");
    printf("Porcentagens do tempo total de processamento \n");
    printf("PCT do tempo de processamento [s]...: %f \n", 100*(soma_total/soma_total));
    // printf("PCT do tempo total de interpolacao [ms]...: %f \n", 100*(soma_interpolacao/soma_total));
    // printf("PCT do tempo total de decimacao [ms]...: %f \n", 100*(soma_decimacao/soma_total));
    // printf("PCT do tempo total de transferencia de memoria H/D [ms]...: %f \n", 100*(soma_thd/soma_total));
    // printf("PCT do tempo total de transferencia de memoria D/H [ms]...: %f \n", 100*(soma_tdh/soma_total));


    //Liberando a memória alocada
    free(amostra_in);
    free(amostra_out);

    return 0;
}