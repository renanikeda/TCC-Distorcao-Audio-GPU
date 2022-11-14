Ressaltamos que os programas contidos nesse arquivo compactado são excluisivos de nossa autoria,
a nao ser onde executadas funções de outras bibliotecas como a libsndfile e CUDA
    - Gustavo Faia Fagundes
    - Leonardo Posso Benetti
    - Renan Ikeda Fernandes

A estrutura desse arquivo é a seguinte:

Entrega
|-> AudioGuitarra
    |- TCC.cu: programa que gera a distorcao em um audio de guitarra
|
|-> Matlab: scripts Matlab
    |- SimulaJetsonDistorcao.m: simula o processamento da GPU
    |- VerificaAliasing_mais_filtrado.m: estima o aliasing de um estágio de distorcao
    |- VerificaCos.m: aplica algoritmo de calcular as harmonicas da funcao atan()
    |- VerificaGanho.m: compara diferentes graus de distorcao da funcao atan()
    |- VerificaTom.m: plota o espectro da soma do atan() de dois tons
|
|-> ProgramasCPU
    |-> ComInterpolacao
        |- TCC_cpu.cpp: executa um programa na CPU equivalente ao da GPU desenvolvido no projeto
    |-> SemInterpolacao
        |- TCC_sem_interpolacao_cpu.cpp: executa um programa na CPU que implementa apenas as distorcoes, sem interpolacoes
    |-> TimingCPU
        |- TCC_cpu_delay.cpp: programa para levantamento de dados experimentais para delay entre amostra de entrada e saida no caso da CPU
        |- TCC_cpu_time.cpp: programa para levantamento de dados experimentais para tempo de execução do programa no caso da CPU
|
|-> TimingGPU
    |- TCC_gpu_delay_cuda.cu: programa para levantamento de dados experimentais para delay entre amostra de entrada e saida no caso da GPU
    |- TCC_gpu_time.cu: programa para levantamento de dados experimentais para tempo de execução do programa no caso da GPU
|


Para compilar arquivos .cu:
	nvcc -o <nome_executável> <nome_do_script.cu> -I/usr/local/include -lsndfile

Para compilar arquivos .cpp:
	g++ -o <nome_executável> <nome_do_script.cpp> -I/usr/local/include -lsndfile