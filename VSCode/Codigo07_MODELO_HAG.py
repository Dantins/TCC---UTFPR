import tclab
import numpy as np
import time
import matplotlib.pyplot as plt

# Conexão com o TCLab
a = tclab.TCLab()

# Tempo de execução em minutos
run_time = 20

# Número de ciclos (1 ciclo por segundo)
loops = int(60.0 * run_time)
tm = np.zeros(loops)            # Vetor para armazenar o tempo
T1 = np.zeros(loops)            # Temperatura real
T1_ModeloHAG = np.zeros(loops)    # Temperatura modelo HAG
Q1 = np.ones(loops) * 50        # Potência fixa em 50%

a.Q1(50)    # Aplicar potência ao hardware
Q = 50      # Aplicar potência aos modelos

# Gráfico
plt.figure(figsize=(12, 8))
plt.ion()
plt.show()

# Loop principal
start_time = time.time()
prev_time = start_time
dt = 1.0  # intervalo de tempo em segundos

try:
    for i in range(loops):
        # Tempo de espera
        sleep_time = max(0.01, dt - (time.time() - prev_time))
        time.sleep(sleep_time)

        # Atualizar tempo
        t = time.time()
        prev_time = t
        tm[i] = t - start_time

        # Ler temperatura real do TCLab
        T1[i] = a.T1

        # --- Modelo Hägglund ---
        if i <= 0:
            T1_ModeloHAG[i] = T1[0]  # condição inicial
        else:
            T1_ModeloHAG[i] = T1_ModeloHAG[i-1] + 0.00328 * np.exp(-tm[i - 15] / 180) * Q * np.heaviside(i - 15, 0)


        # Atualizar gráfico
        plt.clf()
        ax1 = plt.subplot(2, 1, 1)
        ax1.grid()
        plt.plot(tm[:i+1], T1[:i+1], 'r-', label='Temperatura Real', linewidth=2)
        plt.plot(tm[:i+1], T1_ModeloHAG[:i+1], 'g--', label='Modelo Hägglund', linewidth=2)
        plt.ylabel('Temperatura (°C)', fontsize=14)
        plt.legend(loc='best', fontsize=12)

        ax2 = plt.subplot(2, 1, 2)
        ax2.grid()
        plt.plot(tm[:i+1], Q1[:i+1], 'r-', label='Potência (%)', linewidth=2)
        plt.ylabel('Potência (%)', fontsize=14)
        plt.xlabel('Tempo (s)', fontsize=14)
        plt.legend(loc='best', fontsize=12)

        plt.tight_layout()
        plt.draw()
        plt.pause(0.05)

    # --- Cálculo do Erro e Salvamento ---
    ##erro_medio_NL = np.mean(np.abs(T1 - T1_ModeloNL))
    erro_medio_HAG = np.mean(np.abs(T1 - T1_ModeloHAG))
    ##print(f'\nErro Médio Absoluto entre Temperatura Real e Modelo Não Linear: {erro_medio_NL:.2f} °C')
    print(f'Erro Médio Absoluto entre Temperatura Real e Modelo Hägglund: {erro_medio_HAG:.2f} °C')

    # Salvar gráfico em PNG
    plt.savefig('Grafico_TCLab_Modelo_HAG.png', dpi=300, bbox_inches='tight')

    # Salvar dados em .txt
    np.savetxt('Dados_simulacao_07_HAG.txt', np.column_stack((tm, T1, T1_ModeloHAG, Q1)),
               header='Tempo (s)\tTemperatura Real (°C)\tModelo Hägglund (°C)\tPotência (%)', fmt='%.2f', delimiter='\t')

except KeyboardInterrupt:
    print('Interrompido pelo usuário.')

except Exception as e:
    print(f'Erro: {e}')
    raise

finally:
    a.Q1(0)    # Desliga o aquecedor
    a.LED(0)   # Apaga LED
    a.close()
    print('Simulação concluída.')
