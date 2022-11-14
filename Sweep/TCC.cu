//Bibliotecas
////////////////////////////////////////////////////////////////////////////////////////////////
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include <device_launch_parameters.h>
#include "sndfile.h"
#include <math.h>
#include <cuda_runtime_api.h> 
#include <time.h> 
////////////////////////////////////////////////////////////////////////////////////////////////

//Defines ajustaveis
////////////////////////////////////////////////////////////////////////////////////////////////
#define L 30
#define N_interpol 6
#define N_decimacao 6
#define N_aliasing 7
#define N_operacoes_case0 125//Numero de filtragens no kernel de Interpolacao
#define N_operacoes_case1 125 //Numero de filtragens no kernel de decimacao
#define N_operacoes_case2  125//Numero de filtragens no kernel de distorcao
#define N_filtragens_interpol 4 //Numero de filtragens pelo filtro de interpolacao
#define N_filtragens_decimacao 4//Numero de filtragens pelo filtro de decimacao
#define N_estagios 16 //Numero de estagios (incluindo Interpolacao e Decimacao)
#define N_buffer 10000
////////////////////////////////////////////////////////////////////////////////////////////////

//Defines em função
////////////////////////////////////////////////////////////////////////////////////////////////
#define N_eventos (N_estagios+1)
#define N_streams (N_estagios+2)
#define N_estagios_menos_1 (N_estagios-1)
#define N_estagios_menos_2 (N_estagios-2)
#define N_aliasing_mais_1 (N_aliasing+1)
#define N_interpol_mais_1 (N_interpol+1)
#define N_decimacao_mais_1 (N_decimacao+1)
#define offset_decimacao (N_operacoes_case0+N_operacoes_case2*(N_estagios_menos_2))
#define N_distorcao (N_operacoes_case0+N_operacoes_case1+N_operacoes_case2*(N_estagios_menos_2))
#define N_eventos_mais_um (N_eventos+1)
/////////////////////////////////////////////////////////////////////////////////////////////////

//Kernel de Interpolacao (N_operacoes_case0+N_filtragens_interpol)
__global__ void Interpolacao(
    double* amostra_in, double* amostra_out,
    double* entrada_interpol, double* saida_interpol,
    double* entrada_aliasing, double* saida_aliasing,
    double *a_interpol, double *b_interpol,
    double *a_aliasing, double *b_aliasing,
    double grau_distorcao)

{    
    int gid = blockIdx.x*blockDim.x+threadIdx.x;
    __shared__ int posicao[L];
    __shared__ double soma_entrada[L];
    __shared__ double soma_saida[L];
    double entrada_multiplicada_interpolacao[N_interpol_mais_1];
    double saida_multiplicada_interpolacao[N_interpol];
    double entrada_multiplicada_aliasing[N_aliasing_mais_1];
    double saida_multiplicada_aliasing[N_aliasing];

    posicao[gid]=-gid;
    amostra_out[gid]=amostra_in[0];

    while(posicao[L-1]<N_operacoes_case0+N_filtragens_interpol)
    {    
        soma_entrada[gid]=0;
        soma_saida[gid]=0;

    
        if(posicao[gid]>-1 && posicao[gid]<N_filtragens_interpol) //Interpolacao
        {
            for(int i=N_interpol; i>0; i--)
            {
                entrada_interpol[i+(N_interpol+1)*posicao[gid]]=entrada_interpol[i-1+(N_interpol+1)*posicao[gid]];
            }
            entrada_interpol[(N_interpol+1)*posicao[gid]]=amostra_out[gid];

            for(int i=0; i<N_interpol_mais_1; i++)
            {
                entrada_multiplicada_interpolacao[i]=entrada_interpol[i+(N_interpol+1)*posicao[gid]]*b_interpol[i];
                soma_entrada[gid]+=entrada_multiplicada_interpolacao[i];
            }
            for(int i=0; i<N_interpol; i++)
            {
                saida_multiplicada_interpolacao[i]=saida_interpol[i+(N_interpol)*posicao[gid]]*a_interpol[i];
                soma_saida[gid]+=saida_multiplicada_interpolacao[i];
            }

            amostra_out[gid]=soma_entrada[gid]-soma_saida[gid];

            for(int i=N_interpol-1; i>0; i--)
            {
                saida_interpol[i+(N_interpol)*posicao[gid]]=saida_interpol[i-1+(N_interpol)*posicao[gid]];
            }
            saida_interpol[(N_interpol)*posicao[gid]]=amostra_out[gid];

            if(posicao[gid]==N_filtragens_interpol-1)
            {
                if(amostra_out[gid]>1)
                {
                    amostra_out[gid]=1;
                }
                else if(amostra_out[gid]<-1)
                {
                    amostra_out[gid]=-1;
                }
            }   
        }
        else if(posicao[gid]>(N_filtragens_interpol-1) && posicao[gid]<N_operacoes_case0+(N_filtragens_interpol)) //Distorcao
        {
            for(int i=N_aliasing; i>0; i--)
            {
                entrada_aliasing[i+(posicao[gid]-N_filtragens_interpol)*(N_aliasing_mais_1)]=entrada_aliasing[i-1+(posicao[gid]-N_filtragens_interpol)*(N_aliasing_mais_1)];
            }
            entrada_aliasing[0+(posicao[gid]-N_filtragens_interpol)*(N_aliasing_mais_1)]=atan(grau_distorcao*amostra_out[gid])/atan(grau_distorcao);

            for(int i=0; i<N_aliasing_mais_1; i++)
            {
                entrada_multiplicada_aliasing[i]=entrada_aliasing[i+(posicao[gid]-N_filtragens_interpol)*(N_aliasing_mais_1)]*b_aliasing[i];
                soma_entrada[gid]+=entrada_multiplicada_aliasing[i];
            }
            for(int i=0; i<N_aliasing; i++)
            {
                saida_multiplicada_aliasing[i]=saida_aliasing[i+(posicao[gid]-N_filtragens_interpol)*N_aliasing]*a_aliasing[i];
                soma_saida[gid]+=saida_multiplicada_aliasing[i];
            }

            amostra_out[gid]=soma_entrada[gid]-soma_saida[gid];

            for(int i=N_aliasing-1; i>0; i--)
            {
                saida_aliasing[i+(posicao[gid]-N_filtragens_interpol)*N_aliasing]=saida_aliasing[i-1+(posicao[gid]-N_filtragens_interpol)*N_aliasing];
            }
            saida_aliasing[(posicao[gid]-N_filtragens_interpol)*N_aliasing]=amostra_out[gid];

        }

        posicao[gid]+=1;

        __syncthreads();        
        
            
    }

}

//Kernel de Distorcao (N_operacoes_case2)
__global__ void Distorcao(
    double* amostra_in, double* amostra_out,
    double* entrada_aliasing, double* saida_aliasing,
    double *a_aliasing, double *b_aliasing,
    double grau_distorcao, int estagio
)
{    
    int gid = blockIdx.x*blockDim.x+threadIdx.x;
    __shared__ int posicao[L];
    __shared__ double soma_entrada[L];
    __shared__ double soma_saida[L];
    double entrada_multiplicada_aliasing[N_aliasing_mais_1];
    double saida_multiplicada_aliasing[N_aliasing];
   

    posicao[gid]=-gid;

    while(posicao[L-1]<N_operacoes_case2)
    {    
        soma_entrada[gid]=0;
        soma_saida[gid]=0;

      
        if(posicao[gid]>-1 && posicao[gid]<N_operacoes_case2) //Distorcao
        {
            for(int i=N_aliasing; i>0; i--)
            {
                entrada_aliasing[i+(posicao[gid]+N_operacoes_case0+(estagio-1)*N_operacoes_case2)*(N_aliasing_mais_1)]=entrada_aliasing[i-1+(posicao[gid]+N_operacoes_case0+(estagio-1)*N_operacoes_case2)*(N_aliasing_mais_1)];
            }
            entrada_aliasing[0+(posicao[gid]+N_operacoes_case0+(estagio-1)*N_operacoes_case2)*(N_aliasing_mais_1)]=atan(grau_distorcao*amostra_in[gid])/atan(grau_distorcao);

            for(int i=0; i<N_aliasing_mais_1; i++)
            {
                entrada_multiplicada_aliasing[i]=entrada_aliasing[i+(posicao[gid]+N_operacoes_case0+(estagio-1)*N_operacoes_case2)*(N_aliasing_mais_1)]*b_aliasing[i];
                soma_entrada[gid]+=entrada_multiplicada_aliasing[i];
            }
            for(int i=0; i<N_aliasing; i++)
            {
                saida_multiplicada_aliasing[i]=saida_aliasing[i+(posicao[gid]+N_operacoes_case0+(estagio-1)*N_operacoes_case2)*N_aliasing]*a_aliasing[i];
                soma_saida[gid]+=saida_multiplicada_aliasing[i];
            }

            amostra_in[gid]=soma_entrada[gid]-soma_saida[gid];

            for(int i=N_aliasing-1; i>0; i--)
            {
                saida_aliasing[i+(posicao[gid]+N_operacoes_case0+(estagio-1)*N_operacoes_case2)*N_aliasing]=saida_aliasing[i-1+(posicao[gid]+N_operacoes_case0+(estagio-1)*N_operacoes_case2)*N_aliasing];
            }
            saida_aliasing[(posicao[gid]+N_operacoes_case0+(estagio-1)*N_operacoes_case2)*N_aliasing]=amostra_in[gid];
        }
                 
    
        posicao[gid]+=1;

        __syncthreads();
    }
}

//Kernel de Decimacao (N_operacoes_case1-1 distorcoes)
__global__ void Decimacao(
    double* amostra_in, double* amostra_out,
    double* entrada_aliasing, double* saida_aliasing,
    double* entrada_decimacao, double* saida_decimacao,
    double *a_aliasing, double *b_aliasing,
    double *a_decimacao, double *b_decimacao,
    double grau_distorcao
)
{    
    int gid = blockIdx.x*blockDim.x+threadIdx.x;
    __shared__ int posicao[L];
    __shared__ double soma_entrada[L];
    __shared__ double soma_saida[L];
    double entrada_multiplicada_decimacao[N_decimacao_mais_1];
    double saida_multiplicada_decimacao[N_decimacao];
    double entrada_multiplicada_aliasing[N_aliasing_mais_1];
    double saida_multiplicada_aliasing[N_aliasing];

    posicao[gid]=-gid;

    while(posicao[L-1]<N_operacoes_case1+(N_filtragens_decimacao))
    {    
        soma_entrada[gid]=0;
        soma_saida[gid]=0;

        if(posicao[gid]>N_operacoes_case1-1 && posicao[gid]<N_operacoes_case1+(N_filtragens_decimacao)) //Decimacao
        {
            for(int i=N_decimacao; i>0; i--)
            {
                entrada_decimacao[i+(posicao[gid]-N_operacoes_case1)*(N_decimacao_mais_1)]=entrada_decimacao[i-1+(posicao[gid]-N_operacoes_case1)*(N_decimacao_mais_1)];
            }
            entrada_decimacao[(posicao[gid]-N_operacoes_case1)*(N_decimacao_mais_1)]=amostra_in[gid];

            for(int i=0; i<N_decimacao_mais_1; i++)
            {
                entrada_multiplicada_decimacao[i]=entrada_decimacao[i+(posicao[gid]-N_operacoes_case1)*(N_decimacao_mais_1)]*b_decimacao[i];
                soma_entrada[gid]+=entrada_multiplicada_decimacao[i];
            }
            for(int i=0; i<N_decimacao; i++)
            {
                saida_multiplicada_decimacao[i]=saida_decimacao[i+(posicao[gid]-N_operacoes_case1)*(N_decimacao)]*a_decimacao[i];
                soma_saida[gid]+=saida_multiplicada_decimacao[i];
            }

            amostra_in[gid]=soma_entrada[gid]-soma_saida[gid];

            for(int i=N_decimacao-1; i>0; i--)
            {
                saida_decimacao[i+(posicao[gid]-N_operacoes_case1)*(N_decimacao)]=saida_decimacao[i-1+(posicao[gid]-N_operacoes_case1)*(N_decimacao)];
            }
            saida_decimacao[(posicao[gid]-N_operacoes_case1)*(N_decimacao)]=amostra_in[gid];

            if(gid==0 && posicao[gid]==N_operacoes_case1+(N_filtragens_decimacao)-1)
            {
                amostra_out[0]=amostra_in[gid];
            }

        }
        else if(posicao[gid]>-1 && posicao[gid]<N_operacoes_case1) //Distorcao
        {
            for(int i=N_aliasing; i>0; i--)
            {
                entrada_aliasing[i+(posicao[gid]+offset_decimacao)*(N_aliasing_mais_1)]=entrada_aliasing[i-1+(posicao[gid]+offset_decimacao)*(N_aliasing_mais_1)];
            }
            entrada_aliasing[(posicao[gid]+offset_decimacao)*(N_aliasing_mais_1)]=atan(grau_distorcao*amostra_in[gid])/atan(grau_distorcao);

            for(int i=0; i<N_aliasing_mais_1; i++)
            {
                entrada_multiplicada_aliasing[i]=entrada_aliasing[i+(posicao[gid]+offset_decimacao)*(N_aliasing_mais_1)]*b_aliasing[i];
                soma_entrada[gid]+=entrada_multiplicada_aliasing[i];
            }
            for(int i=0; i<N_aliasing; i++)
            {
                saida_multiplicada_aliasing[i]=saida_aliasing[i+(posicao[gid]+offset_decimacao)*N_aliasing]*a_aliasing[i];
                soma_saida[gid]+=saida_multiplicada_aliasing[i];
            }

            amostra_in[gid]=soma_entrada[gid]-soma_saida[gid];

            for(int i=N_aliasing-1; i>0; i--)
            {
                saida_aliasing[i+(posicao[gid]+offset_decimacao)*N_aliasing]=saida_aliasing[i-1+(posicao[gid]+offset_decimacao)*N_aliasing];
            }
            saida_aliasing[(posicao[gid]+offset_decimacao)*N_aliasing]=amostra_in[gid];
        }

       
        posicao[gid]+=1;

        __syncthreads();
    }
    
}

//Kernel para verificar uma variavel do device 
__global__ void Verifica_buffer(double *bufferd, int i, int nome)
{
    if(nome == 0)
        printf("H/D buffer[%d] = %e\n", i, bufferd[i]);
    else
        printf("D/H buffer[%d] = %e\n", i, bufferd[i]);
}

//Kernel para verificar um vetor do device
__global__ void Verifica_vetor(double *bufferd, int N)
{
    for (int i=0; i<N;i++)
    {
        printf("vetor[%d]=%e\n", i, bufferd[i]);
    }
}

int main(){

    printf("N_distorcao: %d\n", N_distorcao);
    printf("Definindo arquivo de entrada e saída: \n");
    
    //Arquivo de entrada e saida
    ////////////////////////////////////////////////////////////////////////////////////////////////
    //Definições do arquivo de leitura e saida
    SNDFILE *file_in,*file_out ;   //Arquivo de entrada e saída
    SF_INFO sfinfo_in,sfinfo_out ; //Arquivo de informações de entrada e saída
    sfinfo_in.format = 0;          //Documentação do libsnd manda fazer isso para arquivos de leitura
    file_in = sf_open ("sweep_48khz_5s.wav", SFM_READ, &sfinfo_in); //Determina o arquivo de entrada
    sfinfo_out=sfinfo_in;
    file_out = sf_open ("sweep_48khz_5s_2000_0.045.wav", SFM_WRITE, &sfinfo_out); //Determina o arquivo de saída
    sf_command (file_out, SFC_SET_CLIPPING, NULL, SF_TRUE) ;
    printf("Informações sobre o arquivo de entrada: \n");
    printf("Taxa de amostragem = %d , Frames = % d , Canais = % d \n" , (int) sfinfo_in.samplerate,(int) sfinfo_in.frames,(int)sfinfo_in.channels); 
    //Mostra algumas caratersiticas do arquivo de entrada
    int read_count = 1;
    ////////////////////////////////////////////////////////////////////////////////////////////////

    printf("Alocando memória: \n");

    //Definicao das variaveis e alocacao de memoria
    ////////////////////////////////////////////////////////////////////////////////////////////////
    //Definições do Projeto
    const double grau_distorcao=0.045; //Grau de distorcao do atan
    const int N_amostras= sfinfo_in.frames; //Numero de amostras do arquivo de entrada
    //const int N_amostras = 1000; //Para verficacao da timeline

    //Variaveis da CPU
    double *amostra_in_host; //Amostra para a entrada no processamento feito na GPU
    double *amostra_out_host; //Amostra para a saída no processamento feito na GPU

    //Alocação das Variaveis da CPU
    cudaHostAlloc((void**)&amostra_in_host, sizeof(double)  *N_buffer,cudaHostAllocDefault);
    cudaHostAlloc((void**)&amostra_out_host, sizeof(double)  *N_buffer,cudaHostAllocDefault);

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

    //Variaveis da GPU
    double *buffer_device_in; //Buffer do devide para leitura do arquivo
    double *buffer_device_out; //Buffer do para escrever amostras processadas
    double *amostra_estagio; //Buffer intermediario das distorcoes
    double *a_interpol; //Vetor dos coeficientes recursivos do filtro de interpolação 
    double *b_interpol; //Vetor dos coeficientes relacionados a entrada do filtro de interpolação
    double *a_decimacao; //Vetor dos coeficientes recursivos do filtro de decimacao
    double *b_decimacao; //Vetor dos coeficientes relacionados a entrada do filtro de decimação
    double *a_aliasing; //Vetor dos coeficientes recursivos do filtro pós distorção 
    double *b_aliasing; //Vetor dos coeficientes relacionados a entrada do filtro pós distorção
    double *entrada_interpol; //Vetor de entrada na interpolação
    double *saida_interpol; //Vetor de saída na interpolação
    double *entrada_decimacao; //Vetor de entrada na decimação
    double *saida_decimacao; //Vetor de saída na decimação
    double *entrada_aliasing; //Vetor de entrada no filtro pós distorção
    double *saida_aliasing; //Vetor de saída no filtro pós distorção


    //Alocação das Variaveis da GPU
    cudaMalloc((void**)&amostra_estagio, sizeof(double)* L*N_buffer);
    cudaMalloc((void**)&buffer_device_in, sizeof(double)*N_buffer);
    cudaMalloc((void**)&buffer_device_out, sizeof(double)*N_buffer);
    cudaMalloc((void**)&a_interpol, sizeof(double) * N_interpol);
    cudaMalloc((void**)&b_interpol, sizeof(double) * (N_interpol_mais_1));
    cudaMalloc((void**)&a_decimacao, sizeof(double) * N_decimacao);
    cudaMalloc((void**)&b_decimacao, sizeof(double) * (N_decimacao_mais_1));
    cudaMalloc((void**)&a_aliasing, sizeof(double) * N_aliasing);
    cudaMalloc((void**)&b_aliasing, sizeof(double) * (N_aliasing_mais_1));
    cudaMalloc((void**)&saida_interpol, sizeof(double) * (N_filtragens_interpol)*N_interpol);
    cudaMalloc((void**)&entrada_interpol, sizeof(double) * (N_filtragens_interpol)*(N_interpol_mais_1));
    cudaMalloc((void**)&saida_decimacao, sizeof(double) * (N_filtragens_decimacao)*N_decimacao);
    cudaMalloc((void**)&entrada_decimacao, sizeof(double) * (N_filtragens_decimacao)*(N_decimacao_mais_1));
    cudaMalloc((void**)&entrada_aliasing, sizeof(double) * N_distorcao*(N_aliasing_mais_1));
    cudaMalloc((void**)&saida_aliasing, sizeof(double) * N_distorcao*N_aliasing);
    

    //Definição dos coeficientes dos filtros na GPU;
    cudaMemcpy(a_interpol, a_interpol_host, sizeof(double)*N_interpol, cudaMemcpyHostToDevice);
    cudaMemcpy(b_interpol, b_interpol_host, sizeof(double)*(N_interpol_mais_1), cudaMemcpyHostToDevice);
    cudaMemcpy(a_decimacao, a_decimacao_host, sizeof(double)*N_decimacao, cudaMemcpyHostToDevice);
    cudaMemcpy(b_decimacao, b_decimacao_host, sizeof(double)*(N_decimacao_mais_1), cudaMemcpyHostToDevice);
    cudaMemcpy(a_aliasing, a_aliasing_host, sizeof(double)*N_aliasing, cudaMemcpyHostToDevice);
    cudaMemcpy(b_aliasing, b_aliasing_host, sizeof(double)*(N_aliasing_mais_1), cudaMemcpyHostToDevice);


    //Definindo filtros com condicoes iniciais nulas
    cudaMemset((void**)&saida_interpol, 0, sizeof(double)*(N_filtragens_interpol)*N_interpol); 
    cudaMemset((void**)&entrada_interpol, 0, sizeof(double)*(N_filtragens_interpol)*(N_interpol_mais_1)); 
    cudaMemset((void**)&saida_decimacao, 0, sizeof(double)*(N_filtragens_decimacao)*N_decimacao);
    cudaMemset((void**)&entrada_decimacao, 0, sizeof(double)*(N_filtragens_decimacao)*(N_decimacao_mais_1));
    cudaMemset((void**)&saida_aliasing, 0, sizeof(double)*N_aliasing*N_distorcao);
    cudaMemset((void**)&entrada_aliasing, 0, sizeof(double)*(N_aliasing_mais_1)*N_distorcao);
    cudaMemset((void**)&buffer_device_out, 2, sizeof(double)*N_buffer);
    cudaMemset((void**)&buffer_device_in, 0, sizeof(double)*N_buffer);
    
    //Definindo Streams
    cudaStream_t estagio[N_streams];
    for(int i=0; i<(N_streams); i++)
    {
        cudaStreamCreate(&estagio[i]);
    }

    //Definindo Eventos
    cudaEvent_t *fim_estagio;
    fim_estagio=(cudaEvent_t*)malloc(sizeof(cudaEvent_t)* (N_buffer+1)*(N_eventos_mais_um));
    for(int i=0; i<(N_buffer+1)*(N_eventos_mais_um); i++)
    {
        cudaEventCreateWithFlags(&fim_estagio[i],cudaEventDisableTiming);
    }

    //Variavel de confirmacao de eventos (mem dinamica)
    // int* event_ok;
    // event_ok=(int*)calloc((N_buffer+1)*(N_eventos_mais_um),sizeof(int));

    //Variavel de confirmacao de eventos (mem estatica)
    int event_ok[(N_buffer+1)*(N_eventos_mais_um)];
    for(int i=0; i<(N_buffer+1)*(N_eventos_mais_um) ; i++)
    {
        event_ok[i]=0;
    }

    //Variaveis de medida de tempo
    clock_t inicio_total;
	clock_t fim_total;
    double soma_total = 0;
    int pct=0;
    int amostras_enviadas=0;

    //Contadores de cada estagio
    int amostras_processadas[N_streams];
    int amostras_processadas_totais=0;
    for(int i=0; i<N_streams; i++)
    {
        amostras_processadas[i]=0;
    }
    
    //Arquivo de log
    FILE *fp = fopen("log_TCC_jerson_query_debug.txt","w");
 
    ////////////////////////////////////////////////////////////////////////////////////////////////

    printf("Começando processamento: \n");

    inicio_total = clock();
    while (amostras_processadas_totais< N_amostras)
    {
        if(amostras_processadas[0]< N_buffer && amostras_processadas[0]==amostras_processadas[1])
        {
            read_count = (int) sf_read_double (file_in, &amostra_in_host[amostras_processadas[0]], 1);
            cudaStreamSynchronize(estagio[1]);
            if(read_count==1)
            {
                cudaMemcpyAsync(&buffer_device_in[amostras_processadas[0]], &amostra_in_host[amostras_processadas[0]], sizeof(double), cudaMemcpyHostToDevice,estagio[0]);
                cudaEventRecord(fim_estagio[(N_eventos_mais_um)*amostras_processadas[0]],estagio[0]);
                event_ok[(N_eventos_mais_um)*amostras_processadas[0]]=1;
                amostras_processadas[0]++;
                amostras_enviadas++;
            }
        }

        if(cudaEventQuery(fim_estagio[(amostras_processadas[1])*(N_eventos_mais_um)])==cudaSuccess && (event_ok[(amostras_processadas[1])*(N_eventos_mais_um)]==1))
        {
            event_ok[(amostras_processadas[1])*(N_eventos_mais_um)]=0;

            Interpolacao<<<1,L,0,estagio[1]>>>(
                &buffer_device_in[amostras_processadas[1]],&amostra_estagio[L*amostras_processadas[1]],
                entrada_interpol,saida_interpol,
                entrada_aliasing,saida_aliasing,
                a_interpol, b_interpol,
                a_aliasing, b_aliasing,
                grau_distorcao);
            cudaEventRecord(fim_estagio[1+amostras_processadas[1]*(N_eventos_mais_um)],estagio[1]);
            // Verifica_vetor<<<1,1,0,estagio[1]>>>(&amostra_estagio[L*amostras_processadas[1]],L);
            // cudaStreamSynchronize(estagio[1]);
            event_ok[(N_eventos_mais_um)*amostras_processadas[1]+1]=1;
            amostras_processadas[1]++;
            fprintf(fp,"Interpolacao amostras_processadas[1]=%d \n", amostras_processadas[1]);
        }

        for(int z=1; z<N_estagios_menos_1;z++)
        {
            if(cudaEventQuery(fim_estagio[z+amostras_processadas[(z+1)]*(N_eventos_mais_um)])==cudaSuccess && event_ok[z+amostras_processadas[(z+1)]*(N_eventos_mais_um)]==1)
            {
                event_ok[z+amostras_processadas[(z+1)]*(N_eventos_mais_um)]=0;

                Distorcao<<<1,L,0,estagio[(z+1)]>>>(
                    &amostra_estagio[L*amostras_processadas[(z+1)]],&amostra_estagio[L*amostras_processadas[(z+1)]],
                    entrada_aliasing,saida_aliasing,
                    a_aliasing, b_aliasing,
                    grau_distorcao,z);
                cudaEventRecord(fim_estagio[z+1+amostras_processadas[(z+1)]*(N_eventos_mais_um)],estagio[(z+1)]);
                event_ok[z+1+(N_eventos_mais_um)*amostras_processadas[(z+1)]]=1;
                amostras_processadas[(z+1)]++;
                fprintf(fp,"Distorcao amostras_processadas[2]=%d \n", amostras_processadas[(z+1)]);
            }
        } 
       
      
        if(cudaEventQuery(fim_estagio[(N_estagios_menos_1)+amostras_processadas[(N_estagios)]*(N_eventos_mais_um)])==cudaSuccess && event_ok[(N_estagios_menos_1)+amostras_processadas[(N_estagios)]*(N_eventos_mais_um)]==1)
        {
            event_ok[(N_estagios_menos_1)+amostras_processadas[(N_estagios)]*(N_eventos_mais_um)]=0;

            Decimacao<<<1,L,0,estagio[(N_estagios)]>>>(
                &amostra_estagio[L*amostras_processadas[(N_estagios)]],&buffer_device_out[amostras_processadas[(N_estagios)]],
                entrada_aliasing,saida_aliasing,
                entrada_decimacao, saida_decimacao,
                a_aliasing, b_aliasing,
                a_decimacao, b_decimacao,
                grau_distorcao);
            cudaEventRecord(fim_estagio[(N_estagios)+amostras_processadas[(N_estagios)]*(N_eventos_mais_um)],estagio[(N_estagios)]);
            event_ok[N_estagios+amostras_processadas[(N_estagios)]*(N_eventos_mais_um)]=1;
            amostras_processadas[(N_estagios)]++;
            fprintf(fp,"Decimacao amostras_processadas[(N_estagios)]=%d \n", amostras_processadas[(N_estagios)]);
        }
        
        if(cudaEventQuery(fim_estagio[(N_estagios)+(N_eventos_mais_um)*amostras_processadas[(N_eventos)]])==cudaSuccess && event_ok[(N_estagios)+(N_eventos_mais_um)*amostras_processadas[(N_eventos)]]==1)
        {
            event_ok[(N_estagios)+(N_eventos_mais_um)*amostras_processadas[(N_eventos)]]=0;

            if(amostras_processadas_totais % ((int)sfinfo_in.frames / 10) == 0)
            {
                printf("Porcentagem do processamento = %d%c \n", pct,37);
                pct+=10;
            }
            cudaMemcpyAsync(&amostra_out_host[amostras_processadas[(N_eventos)]], &buffer_device_out[amostras_processadas[(N_eventos)]], sizeof(double), cudaMemcpyDeviceToHost,estagio[(N_eventos)]);
            cudaEventRecord(fim_estagio[(N_eventos)+(N_eventos_mais_um)*amostras_processadas[(N_eventos)]],estagio[(N_eventos)]);
            event_ok[(N_eventos)+(N_eventos_mais_um)*amostras_processadas[(N_eventos)]]=1;
            amostras_processadas[(N_eventos)]++;
            amostras_processadas_totais++;
        }


        if(amostras_processadas[(N_eventos)] == (N_buffer) && amostras_processadas[0]== (N_buffer) )
        {
            if(cudaEventQuery(fim_estagio[(N_eventos)+(N_eventos_mais_um)*(amostras_processadas[(N_eventos)]-1)])==cudaSuccess && event_ok[(N_eventos)+(N_eventos_mais_um)*(amostras_processadas[(N_eventos)]-1)]==1)
            {
                sf_write_double (file_out, amostra_out_host, N_buffer ) ;
                event_ok[(N_eventos)+(N_eventos_mais_um)*(amostras_processadas[(N_eventos)]-1)]=0;
                for(int k=0; k<N_eventos_mais_um;k++)
                {
                    amostras_processadas[k]=0;
                }
            }
        }

    }   

    cudaDeviceSynchronize();
    
    ////////////////////////////////////////////////////////////////////////////////////////////////

    printf("Fim do processamento: \n");

    fim_total = clock();

    soma_total += ((double)(fim_total - inicio_total)) / CLOCKS_PER_SEC;

    sf_write_double (file_out, amostra_out_host, amostras_processadas[(N_eventos)] ) ;

    //Prints finais 
    ////////////////////////////////////////////////////////////////////////////////////////////////
    printf("Nro de blocos processados.....: %d\n", amostras_processadas_totais);
    printf("\n");
    printf("Tempo Total \n"); 
    printf("Tempo total de processamento com tranf de memoria [s]...: %f \n", soma_total);

    printf("\n");
    printf("Porcentagens do tempo total de processamento \n");
    printf("PCT do tempo de processamento [s]...: %f \n", 100*(soma_total/soma_total));
  

    printf("Liberando a memória alocada: \n");

   
    //Liberando a memória alocada
    ////////////////////////////////////////////////////////////////////////////////////////////////
    // free(event_ok);

    
    cudaFree(amostra_in_host);
    cudaFree(amostra_out_host);

    cudaFree(buffer_device_in);
    cudaFree(buffer_device_out);
    cudaFree(amostra_estagio);

    cudaFree(a_interpol);
    cudaFree(b_interpol);
    cudaFree(a_decimacao);
    cudaFree(b_decimacao);
    cudaFree(a_aliasing);
    cudaFree(b_aliasing);
    cudaFree(entrada_interpol);
    cudaFree(saida_interpol);
    cudaFree(entrada_decimacao);
    cudaFree(saida_decimacao);
    cudaFree(entrada_aliasing);
    cudaFree(saida_aliasing);

    for(int i=0; i<(N_buffer+1)*(N_eventos); i++)
    {
        cudaEventDestroy(fim_estagio[i]);
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////

     printf("Fim: \n");
 }