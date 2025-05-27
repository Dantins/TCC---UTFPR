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
T1_ModeloNL = np.zeros(loops)   # Temperatura do modelo não linear
Q1 = np.ones(loops) * 50        # Potência fixa em 50%

# Aplicar potência inicial
a.Q1(50)

# Criar gráfico
plt.figure(figsize=(12, 8))
plt.ion()
plt.show()

# Loop principal
start_time = time.time()
prev_time = start_time

try:
    for i in range(loops):
        # Tempo de espera
        sleep_time = max(0.01, 1.0 - (time.time() - prev_time))
        time.sleep(sleep_time)

        # Atualizar tempo
        t = time.time()
        prev_time = t
        tm[i] = t - start_time

        # Ler temperatura real do TCLab
        T1[i] = a.T1

        # Informações para o modelo
        Taquecedor = T1_ModeloNL[i-1] + 273.15 if i > 0 else a.T1 + 273.15  # Temperatura simulada em Kelvin
        Tambiente = 301.15  # Temperatura ambiente (28°C em Kelvin)
        Alpha = 0.01
        Cp = 500
        A = 0.0012
        m = 0.004
        U = 10
        Emissividade = 0.9
        Boltzmann = 5.67e-8
        Q = 50
        dt = 1.0  # intervalo de tempo em segundos

        # Modelo não linear: cálculo de dT/dt
        dTdt = ((Alpha / (m * Cp)) * Q) + ((U * A) / (m * Cp)) * (Tambiente - Taquecedor) + ((Emissividade * Boltzmann * A) / (m * Cp)) * (Tambiente**4 - Taquecedor**4)

        if i == 0:
            T1_ModeloNL[i] = a.T1  # Condição inicial: T real no instante 0
        else:
            T1_ModeloNL[i] = T1_ModeloNL[i-1] + dTdt * dt  # Temperatura estimada

        # Atualizar gráfico
        plt.clf()
        ax1 = plt.subplot(2, 1, 1)
        ax1.grid()
        plt.plot(tm[:i+1], T1[:i+1], 'r-', label='Temperatura Real', linewidth=2)
        plt.plot(tm[:i+1], T1_ModeloNL[:i+1], 'b--', label='Modelo Não Linear', linewidth=2)
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

    # Erro Médio Absoluto
    erro_medio = np.mean(np.abs(T1 - T1_ModeloNL))
    print(f'\nErro Médio Absoluto entre Temperatura Real e Modelo: {erro_medio:.2f} °C')

    # Salvar gráfico em PNG
    plt.savefig('Grafico_TCLab_ModeloNL.png', dpi=300, bbox_inches='tight')

    # Salvar dados em .txt
    np.savetxt('Dados_simulacao_01.txt', np.column_stack((tm, T1, T1_ModeloNL, Q1)),
               header='Tempo (s)\tTemperatura Real (°C)\tTemperatura Modelo (°C)\tPotência (%)', fmt='%.2f', delimiter='\t')

# Capturar interrupção manual (Ctrl+C)
except KeyboardInterrupt:
    print('Interrompido pelo usuário.')

# Capturar outros erros
except Exception as e:
    print(f'Erro: {e}')
    raise

# Finalizar e garantir que o aquecedor desligue
finally:
    a.Q1(0)    # Desliga o aquecedor
    a.LED(0)   # Apaga LED
    a.close()
    print('Simulação concluída.')
