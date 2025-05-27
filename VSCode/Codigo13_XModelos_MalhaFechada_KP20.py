import tclab
import numpy as np
import time
import matplotlib.pyplot as plt

# Conexão com o TCLab
a = tclab.TCLab()

# Tempo de execução em minutos
run_time = 20
loops = int(60.0 * run_time)

# Inicialização dos vetores
tm = np.zeros(loops)                #Tempo
T1 = np.zeros(loops)                #Temperatura Real
T1_ModeloNL = np.zeros(loops)
T1_ModeloL = np.zeros(loops)
T1_ModeloHAG = np.zeros(loops)
T1_ModeloSmith = np.zeros(loops)
T1_ModeloSundaresan = np.zeros(loops)
T1_ModeloZN = np.zeros(loops)
Q1 = np.zeros(loops)                #Potência Real Aplicada
Q_NL = np.zeros(loops)              #Potência Modelo Não Linear
Q_L = np.zeros(loops)               #Potência Modelo Linear
Q_HAG = np.zeros(loops)             #Potência Modelo HAG
Q_Smith = np.zeros(loops)           #Potência Modelo Smith
Q_Sundaresan = np.zeros(loops)      #Potência Modelo Sundaresan
Q_ZN = np.zeros(loops)              #Potência Modelo ZN

# Setpoint e parâmetros PID
SP = 50  # °C
Kp, Ki, Kd = 20, 0, 0
integral = 0.0
erro_anterior = 0.0

# Parâmetros físicos para o modelo não linear
Alpha = 0.01
Cp = 500
A = 0.0012
m = 0.004
U = 8
Emissividade = 0.9
Boltzmann = 5.67e-8
Tambiente = 301.15  # Kelvin
dt = 1              # segundo
L = 15              # Atraso de transporte

# Criar gráfico interativo
plt.figure(figsize=(10, 7))
plt.ion()
plt.show()

start_time = time.time()
prev_time = start_time

try:
    for i in range(loops):
        # Controle de tempo
        sleep_time = max(0.01, 1.0 - (time.time() - prev_time))
        time.sleep(sleep_time)

        t = time.time()
        prev_time = t
        tm[i] = t - start_time

        # Medição de temperatura
        T1[i] = a.T1

        # Erros
        Erro_real = SP - T1[i]
        Erro_NL = SP - T1_ModeloNL[i-1] if i > 0 else Erro_real
        Erro_L = SP - T1_ModeloL[i-1] if i > 0 else Erro_real
        Erro_HAG = SP - T1_ModeloHAG[i-1] if i > 0 else Erro_real
        Erro_Smith = SP - T1_ModeloSmith[i-1] if i > 0 else Erro_real
        Erro_Sundaresan = SP - T1_ModeloSundaresan[i-1] if i > 0 else Erro_real
        Erro_ZN = SP - T1_ModeloZN[i-1] if i > 0 else Erro_real

        # Controle PID sistema real
        integral += Erro_real * dt
        derivada = (Erro_real - erro_anterior) / dt
        Q_real = Kp * Erro_real + Ki * integral + Kd * derivada
        erro_anterior = Erro_real
        Q1[i] = np.clip(Q_real, 0, 100)  # Limita entre 0 e 100
        a.Q1(Q1[i])                      # Aplica a potência no transistor

        # Controle PID modelos
        Q_NL[i] = np.clip(Kp * Erro_NL, 0, 100)
        Q_L[i] = np.clip(Kp * Erro_L, 0, 100)
        Q_HAG[i] = np.clip(Kp * Erro_HAG, 0, 100)
        Q_Smith[i] = np.clip(Kp * Erro_Smith, 0, 100)
        Q_Sundaresan[i] = np.clip(Kp * Erro_Sundaresan, 0, 100)
        Q_ZN[i] = np.clip(Kp * Erro_ZN, 0, 100)

        # --- Modelo Não Linear ---
        if i < 15:
            T1_ModeloNL[i] = T1[0]
        else:
            Taquecedor = T1_ModeloNL[i-1] + 273.15
            dTdt = ((Alpha / (m * Cp)) * Q_NL[i]) + ((U * A) / (m * Cp)) * (Tambiente - Taquecedor) \
                   + ((Emissividade * Boltzmann * A) / (m * Cp)) * (Tambiente**4 - Taquecedor**4)
            T1_ModeloNL[i] = T1_ModeloNL[i-1] + dTdt * dt

        # Modelo Linear
        if i <= 0:
            T1_ModeloL[i] = T1[0]  # Condição inicial
        else:
            T1_ModeloL[i] = T1_ModeloL[i-1] = T1_ModeloL[i-1] + 0.00499 * np.exp(-tm[i - 15] / 123) * Q_L[i] * np.heaviside(i - 15, 0)

        # --- Modelo ZN ---
        if i <= 0:
            T1_ModeloZN[i] = T1[0]  # condição inicial
        else:
            T1_ModeloZN[i] = T1_ModeloZN[i-1] + 0.00272 * np.exp(-tm[i - 15] / 217) * Q_ZN[i] * np.heaviside(i - 15, 0)

        # --- Modelo Hägglund ---
        if i <= 0:
            T1_ModeloHAG[i] = T1[0]  # condição inicial
        else:
            T1_ModeloHAG[i] = T1_ModeloHAG[i-1] + 0.00328 * np.exp(-tm[i - 15] / 180) * Q_HAG[i] * np.heaviside(i - 15, 0)    

        # --- Modelo Smith ---
        if i <= 0:
            T1_ModeloSmith[i] = T1[0]  # condição inicial
        else:
            T1_ModeloSmith[i] = T1_ModeloSmith[i-1] + 0.00340 * np.exp(-tm[i - 25] / 174) * Q_Smith[i] * np.heaviside(i - 25, 0)    

        # --- Modelo Sundaresan ---
        if i <= 0:
            T1_ModeloSundaresan[i] = T1[0]  # condição inicial
        else:
            T1_ModeloSundaresan[i] = T1_ModeloSundaresan[i-1] + 0.00297 * np.exp(-tm[i - 16] / 199) * Q_Sundaresan[i] * np.heaviside(i - 16, 0)


        # Gráficos
        plt.clf()
        ax1 = plt.subplot(2, 1, 1)
        ax1.grid()
        ax1.plot(tm[:i+1], T1[:i+1], 'r-', label='Temperatura Real')
        ax1.plot(tm[:i+1], T1_ModeloNL[:i+1], 'b--', label='Modelo Não Linear')
        ax1.plot(tm[:i+1], T1_ModeloL[:i+1], 'g-.', label='Modelo Linear')
        ax1.plot(tm[:i+1], T1_ModeloZN[:i+1], 'k-.', label='Modelo ZN')
        ax1.plot(tm[:i+1], T1_ModeloHAG[:i+1], 'c-.', label='Modelo HAG')
        ax1.plot(tm[:i+1], T1_ModeloSmith[:i+1], 'm-.', label='Modelo Smith')
        ax1.plot(tm[:i+1], T1_ModeloSundaresan[:i+1], 'y-.', label='Modelo Sundaresan')
        ax1.axhline(SP, color='gray', linestyle='--', label='Setpoint')
        ax1.set_ylabel('Temperatura (°C)')
        ax1.legend(loc='best')

        ax2 = plt.subplot(2, 1, 2)
        ax2.grid()
        ax2.plot(tm[:i+1], Q1[:i+1], 'r-', label='Potência Aplicada')
        ax2.plot(tm[:i+1], Q_NL[:i+1], 'b--', label='Modelo Não Linear')
        ax2.plot(tm[:i+1], Q_L[:i+1], 'g-.', label='Modelo Linear')
        ax2.plot(tm[:i+1], Q_ZN[:i+1], 'k-.', label='Modelo ZN')
        ax2.plot(tm[:i+1], Q_HAG[:i+1], 'c-.', label='Modelo HAG')
        ax2.plot(tm[:i+1], Q_Smith[:i+1], 'm-.', label='Modelo Smith')
        ax2.plot(tm[:i+1], Q_Sundaresan[:i+1], 'y-.', label='Modelo Sundaresan')
        ax2.set_ylabel('Potência (%)')
        ax2.set_xlabel('Tempo (s)')
        ax2.legend(loc='best')

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
    print(f'Erro Médio Absoluto entre Temperatura Real e Modelo Hägglund: {erro_medio_HAG:.2f} °C')
    print(f'Erro Médio Absoluto entre Temperatura Real e Modelo Smith: {erro_medio_Smith:.2f} °C')
    print(f'Erro Médio Absoluto entre Temperatura Real e Modelo Sundaresan: {erro_medio_Sundaresan:.2f} °C')
    print(f'Erro Médio Absoluto entre Temperatura Real e Modelo Ziegler_Nichols: {erro_medio_ZN:.2f} °C')

    # Salvar gráfico em PNG
    plt.savefig('Grafico_TCLab_Todos_Modelos_Comparados_Malha_Fechada_KP20.png', dpi=300, bbox_inches='tight')

    # Salvar dados em .txt
    np.savetxt('Dados_simulacao_13_XModelos_Comparados_Malha_Fechada_KP20.txt', np.column_stack((tm, T1, T1_ModeloNL, T1_ModeloL, T1_ModeloHAG, T1_ModeloSmith, T1_ModeloSundaresan, T1_ModeloZN, Q1)),
               header='Tempo (s)\tTemperatura Real (°C)\tModelo Não Linear (°C)\tModelo Linear (°C)\tModelo Hägglund (°C)\tModelo Smith (°C)\tModelo Sundaresan (°C)\tModelo ZN (°C)\tPotência (%)', fmt='%.2f', delimiter='\t')

except KeyboardInterrupt:
    print('Interrompido pelo usuário.')

except Exception as e:
    print(f'Erro: {e}')
    raise

finally:
    a.Q1(0)
    a.LED(0)
    a.close()
    print('Simulação concluída.')
