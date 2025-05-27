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
tm = np.zeros(loops)
T1 = np.zeros(loops)
T1_ModeloNL = np.zeros(loops)
T1_ModeloL = np.zeros(loops)
T1_ModeloHAG = np.zeros(loops)
T1_ModeloSmith = np.zeros(loops)
T1_ModeloSundaresan = np.zeros(loops)
T1_ModeloZN = np.zeros(loops)
Q1 = np.zeros(loops)
Q_NL = np.zeros(loops)
Q_L = np.zeros(loops)
Q_HAG = np.zeros(loops)
Q_Smith = np.zeros(loops)
Q_Sundaresan = np.zeros(loops)
Q_ZN = np.zeros(loops)

# Setpoint e parâmetros PID ideal
SP = 50
Kp = 20
Ti = 50
Td = 0
integral = 0.0
erro_anterior = 0.0

integral_NL = integral_L = integral_HAG = integral_Smith = integral_Sundaresan = integral_ZN = 0.0
erro_anterior_NL = erro_anterior_L = erro_anterior_HAG = erro_anterior_Smith = erro_anterior_Sundaresan = erro_anterior_ZN = 0.0

Alpha = 0.01
Cp = 500
A = 0.0012
m = 0.004
U = 8
Emissividade = 0.9
Boltzmann = 5.67e-8
Tambiente = 301.15
dt = 1
L = 15

plt.figure(figsize=(10, 7))
plt.ion()
plt.show()

start_time = time.time()
prev_time = start_time

try:
    for i in range(loops):
        sleep_time = max(0.01, 1.0 - (time.time() - prev_time))
        time.sleep(sleep_time)
        t = time.time()
        prev_time = t
        tm[i] = t - start_time
        T1[i] = a.T1

        Erro_real = SP - T1[i]
        Erro_NL = SP - T1_ModeloNL[i-1] if i > 0 else Erro_real
        Erro_L = SP - T1_ModeloL[i-1] if i > 0 else Erro_real
        Erro_HAG = SP - T1_ModeloHAG[i-1] if i > 0 else Erro_real
        Erro_Smith = SP - T1_ModeloSmith[i-1] if i > 0 else Erro_real
        Erro_Sundaresan = SP - T1_ModeloSundaresan[i-1] if i > 0 else Erro_real
        Erro_ZN = SP - T1_ModeloZN[i-1] if i > 0 else Erro_real

        # --- Controle PID Ideal Real ---
        integral += Erro_real * dt
        derivada = (Erro_real - erro_anterior) / dt
        Q_real = Kp * (Erro_real + (1/Ti) * integral + Td * derivada)
        if Q_real > 100 and Erro_real > 0:
            integral -= Erro_real * dt
        elif Q_real < 0 and Erro_real < 0:
            integral -= Erro_real * dt
        Q1[i] = np.clip(Q_real, 0, 100)
        erro_anterior = Erro_real
        a.Q1(Q1[i])

        # PID Ideal para modelos
        def pid_model_ideal(Erro, erro_ant, integral_acum, Kp, Ti, Td):
            integral_acum += Erro * dt
            derivada = (Erro - erro_ant) / dt
            Q = Kp * (Erro + (1/Ti) * integral_acum + Td * derivada)
            if Q > 100 and Erro > 0:
                integral_acum -= Erro * dt
            elif Q < 0 and Erro < 0:
                integral_acum -= Erro * dt
            return np.clip(Q, 0, 100), integral_acum, Erro

        Q_NL[i], integral_NL, erro_anterior_NL = pid_model_ideal(Erro_NL, erro_anterior_NL, integral_NL, Kp, Ti, Td)
        Q_L[i], integral_L, erro_anterior_L = pid_model_ideal(Erro_L, erro_anterior_L, integral_L, Kp, Ti, Td)
        Q_HAG[i], integral_HAG, erro_anterior_HAG = pid_model_ideal(Erro_HAG, erro_anterior_HAG, integral_HAG, Kp, Ti, Td)
        Q_Smith[i], integral_Smith, erro_anterior_Smith = pid_model_ideal(Erro_Smith, erro_anterior_Smith, integral_Smith, Kp, Ti, Td)
        Q_Sundaresan[i], integral_Sundaresan, erro_anterior_Sundaresan = pid_model_ideal(Erro_Sundaresan, erro_anterior_Sundaresan, integral_Sundaresan, Kp, Ti, Td)
        Q_ZN[i], integral_ZN, erro_anterior_ZN = pid_model_ideal(Erro_ZN, erro_anterior_ZN, integral_ZN, Kp, Ti, Td)

        if i < L:
            T1_ModeloNL[i] = T1[0]
        else:
            Taquecedor = T1_ModeloNL[i-1] + 273.15
            dTdt = ((Alpha / (m * Cp)) * Q_NL[i]) + ((U * A) / (m * Cp)) * (Tambiente - Taquecedor) \
                   + ((Emissividade * Boltzmann * A) / (m * Cp)) * (Tambiente**4 - Taquecedor**4)
            T1_ModeloNL[i] = T1_ModeloNL[i-1] + dTdt * dt

        if i <= 0:
            T1_ModeloL[i] = T1_ModeloZN[i] = T1_ModeloHAG[i] = T1_ModeloSmith[i] = T1_ModeloSundaresan[i] = T1[0]
        else:
            T1_ModeloL[i] = T1_ModeloL[i-1] + 0.00499 * np.exp(-tm[i - L] / 123) * Q_L[i] * np.heaviside(i - L, 0)
            T1_ModeloZN[i] = T1_ModeloZN[i-1] + 0.00272 * np.exp(-tm[i - L] / 217) * Q_ZN[i] * np.heaviside(i - L, 0)
            T1_ModeloHAG[i] = T1_ModeloHAG[i-1] + 0.00328 * np.exp(-tm[i - L] / 180) * Q_HAG[i] * np.heaviside(i - L, 0)
            T1_ModeloSmith[i] = T1_ModeloSmith[i-1] + 0.00340 * np.exp(-tm[i - 25] / 174) * Q_Smith[i] * np.heaviside(i - 25, 0)
            T1_ModeloSundaresan[i] = T1_ModeloSundaresan[i-1] + 0.00297 * np.exp(-tm[i - 16] / 199) * Q_Sundaresan[i] * np.heaviside(i - 16, 0)

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

    plt.savefig('Grafico_TCLab_Todos_Modelos_Comparados_Malha_Fechada_PI_50.png', dpi=300, bbox_inches='tight')
    np.savetxt('Dados_simulacao_PI_50.txt',
               np.column_stack((tm, T1, T1_ModeloNL, T1_ModeloL, T1_ModeloHAG, T1_ModeloSmith, T1_ModeloSundaresan, T1_ModeloZN, Q1)),
               header='Tempo (s)\tTemperatura Real (°C)\tModelo Não Linear (°C)\tModelo Linear (°C)\tModelo Hägglund (°C)\tModelo Smith (°C)\tModelo Sundaresan (°C)\tModelo ZN (°C)\tPotência (%)',
               fmt='%.2f', delimiter='\t')

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
