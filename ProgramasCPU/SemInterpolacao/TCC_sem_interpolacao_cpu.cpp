#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "sndfile.h"
#include <math.h>
#include <time.h>

int main(int argc, char *argv[ ])
{
    // argc -> nro de argumentos na chamada da funcao
    // argv -> vetor com os argumentos
    if (argc != 3)
    {
        printf("Erro na chamada da funcao.\n modelo: \n\n ./nome_programa [Nro_distorcoes (int)] [grau_distorcao(double)]\n");
        return 1;
    }
    int Nro_distorcoes;
    double grau_distorcao;
    Nro_distorcoes = atoi(argv[1]);
    grau_distorcao = atof(argv[2]);
    printf("Nro distorcoes: %d\n", Nro_distorcoes);
    printf("Grau distorcao: %f\n", grau_distorcao);

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
    
    int N_amostras;
    N_amostras = (int)sfinfo_in.frames;

    // Alocacao do vetor de amostras
    double *amostras;
    amostras = (double*)malloc(sizeof(double) * N_amostras);

    clock_t inicio_total;
	clock_t fim_total;
    double soma_total;

    printf("Começando o loop...\n");

    inicio_total = clock();
    
    int read_count;
    read_count = (int) sf_read_double (file_in, amostras, (int)sfinfo_in.frames);

    for(int iter_amostra=0; iter_amostra<N_amostras; iter_amostra++)
    {
        for (int iter_distorcao=0; iter_distorcao<Nro_distorcoes; iter_distorcao++)
        {
            amostras[iter_amostra] = atan(grau_distorcao*amostras[iter_amostra])/atan(grau_distorcao);
        }
    }

    sf_write_double (file_out, amostras, N_amostras) ;
    
    fim_total = clock();
    soma_total += ((double)(fim_total - inicio_total)) / CLOCKS_PER_SEC;

    printf("Tempo total: %f \n", soma_total);

    return 0;
}
