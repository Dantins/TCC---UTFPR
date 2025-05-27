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
T1_ModeloNL = np.zeros(loops)   # Temperatura modelo não linear
T1_ModeloL = np.zeros(loops)    # Temperatura modelo linear
T1_ModeloHAG = np.zeros(loops)    # Temperatura modelo hag
T1_ModeloSmith = np.zeros(loops)    # Temperatura modelo smith
T1_ModeloSundaresan = np.zeros(loops)    # Temperatura modelo sundaresan
T1_ModeloZN = np.zeros(loops)    # Temperatura modelo ZN
Q1 = np.ones(loops) * 50        # Potência fixa em 50%

a.Q1(50)    # Aplicar potência ao hardware
Q = 50      # Aplicar potência aos modelos

# Parâmetros do modelo não linear
Tambiente = 301.15  # 29 °C em Kelvin
Alpha = 0.01
Cp = 500
A = 0.0012
m = 0.004
U = 8
Emissividade = 0.9
Boltzmann = 5.67e-8

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

        # --- Modelo Não Linear ---
        if i < 15:
            T1_ModeloNL[i] = T1[0]
        else:
            Taquecedor = T1_ModeloNL[i-1] + 273.15
            dTdt = ((Alpha / (m * Cp)) * Q) + ((U * A) / (m * Cp)) * (Tambiente - Taquecedor) \
                   + ((Emissividade * Boltzmann * A) / (m * Cp)) * (Tambiente**4 - Taquecedor**4)
            T1_ModeloNL[i] = T1_ModeloNL[i-1] + dTdt * dt

        # Modelo Linear
        if i <= 0:
            T1_ModeloL[i] = T1[0]  # Condição inicial
        else:
            T1_ModeloL[i] = T1_ModeloL[i-1] = T1_ModeloL[i-1] + 0.00499 * np.exp(-tm[i - 15] / 123) * Q * np.heaviside(i - 15, 0)

        # --- Modelo ZN ---
        if i <= 0:
            T1_ModeloZN[i] = T1[0]  # condição inicial
        else:
            T1_ModeloZN[i] = T1_ModeloZN[i-1] + 0.00272 * np.exp(-tm[i - 15] / 217) * Q * np.heaviside(i - 15, 0)

        # --- Modelo Hägglund ---
        if i <= 0:
            T1_ModeloHAG[i] = T1[0]  # condição inicial
        else:
            T1_ModeloHAG[i] = T1_ModeloHAG[i-1] + 0.00328 * np.exp(-tm[i - 15] / 180) * Q * np.heaviside(i - 15, 0)    

        # --- Modelo Smith ---
        if i <= 0:
            T1_ModeloSmith[i] = T1[0]  # condição inicial
        else:
            T1_ModeloSmith[i] = T1_ModeloSmith[i-1] + 0.00340 * np.exp(-tm[i - 25] / 174) * Q * np.heaviside(i - 25, 0)    

        # --- Modelo Sundaresan ---
        if i <= 0:
            T1_ModeloSundaresan[i] = T1[0]  # condição inicial
        else:
            T1_ModeloSundaresan[i] = T1_ModeloSundaresan[i-1] + 0.00297 * np.exp(-tm[i - 16] / 199) * Q * np.heaviside(i - 16, 0)

        # Atualizar gráfico
        plt.clf()
        ax1 = plt.subplot(2, 1, 1)
        ax1.grid()
        plt.plot(tm[:i+1], T1[:i+1], 'r-', label='Temperatura Real', linewidth=2)
        plt.plot(tm[:i+1], T1_ModeloNL[:i+1], 'b--', label='Modelo Não Linear', linewidth=2)
        plt.plot(tm[:i+1], T1_ModeloL[:i+1], 'g--', label='Modelo Linear', linewidth=2)
        plt.plot(tm[:i+1], T1_ModeloSmith[:i+1], 'y--', label='Modelo Smith', linewidth=2)
        plt.plot(tm[:i+1], T1_ModeloHAG[:i+1], 'c--', label='Modelo Hägglund', linewidth=2)
        plt.plot(tm[:i+1], T1_ModeloZN[:i+1], 'm--', label='Modelo Ziegler_Nichols', linewidth=2)
        plt.plot(tm[:i+1], T1_ModeloSundaresan[:i+1], 'k--', label='Modelo Sundaresan', linewidth=2)
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
    erro_medio_NL = np.mean(np.abs(T1 - T1_ModeloNL))
    erro_medio_L = np.mean(np.abs(T1 - T1_ModeloL))
    erro_medio_HAG = np.mean(np.abs(T1 - T1_ModeloHAG))
    erro_medio_Smith = np.mean(np.abs(T1 - T1_ModeloSmith))
    erro_medio_Sundaresan = np.mean(np.abs(T1 - T1_ModeloSundaresan))
    erro_medio_ZN = np.mean(np.abs(T1 - T1_ModeloZN))
    print(f'\nErro Médio Absoluto entre Temperatura Real e Modelo Não Linear: {erro_medio_NL:.2f} °C')
    print(f'Erro Médio Absoluto entre Temperatura Real e Modelo Linear: {erro_medio_L:.2f} °C')
    print(f'Erro Médio Absoluto entre Temperatura Real e Modelo Ziegler_Nichols: {erro_medio_ZN:.2f} °C')
    print(f'Erro Médio Absoluto entre Temperatura Real e Modelo Hägglund: {erro_medio_HAG:.2f} °C')
    print(f'Erro Médio Absoluto entre Temperatura Real e Modelo Smith: {erro_medio_Smith:.2f} °C')
    print(f'Erro Médio Absoluto entre Temperatura Real e Modelo Sundaresan: {erro_medio_Sundaresan:.2f} °C')


    # Salvar gráfico em PNG
    plt.savefig('Grafico_TCLab_Todos_Modelos_Comparados_Malha_Aberta.png', dpi=300, bbox_inches='tight')

    # Salvar dados em .txt
    np.savetxt('Dados_simulacao_10_Modelos_Comparados_Malha_Aberta.txt', np.column_stack((tm, T1, T1_ModeloNL, T1_ModeloL, T1_ModeloHAG, T1_ModeloSmith, T1_ModeloSundaresan, T1_ModeloZN, Q1)),
               header='Tempo (s)\tTemperatura Real (°C)\tModelo Não Linear (°C)\tModelo Linear (°C)\tModelo Hägglund (°C)\tModelo Smith (°C)\tModelo Sundaresan (°C)\tModelo ZN (°C)\tPotência (%)', fmt='%.2f', delimiter='\t')

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
