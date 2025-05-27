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

# Setpoint e tolerância
SP = 50  # °C
tolerancia = 0.5

# Parâmetros físicos para o modelo não linear
Alpha = 0.01
Cp = 500
A = 0.0012
m = 0.004
U = 8
Emissividade = 0.9
Boltzmann = 5.67e-8
Tambiente = 301.15  # Kelvin
dt = 1
L = 15

# Criar gráfico interativo
plt.figure(figsize=(10, 7))
plt.ion()
plt.show()

# Função de controle ON/OFF com tolerância
def controle_onoff(erro, Q_anterior):
    if erro > tolerancia:
        return 100
    elif erro < -tolerancia:
        return 0
    else:
        return Q_anterior

start_time = time.time()
prev_time = start_time

try:
    for i in range(loops):
        sleep_time = max(0.01, 1.0 - (time.time() - prev_time))
        time.sleep(sleep_time)

        t = time.time()
        prev_time = t
        tm[i] = t - start_time

        # Leitura da temperatura real
        T1[i] = a.T1

        # Controle ON/OFF real
        erro_real = SP - T1[i]
        Q1[i] = controle_onoff(erro_real, Q1[i-1] if i > 0 else 0)
        a.Q1(Q1[i])

        # ON/OFF para modelos
        erro_NL = SP - T1_ModeloNL[i-1] if i > 0 else erro_real
        Q_NL[i] = controle_onoff(erro_NL, Q_NL[i-1] if i > 0 else 0)

        erro_L = SP - T1_ModeloL[i-1] if i > 0 else erro_real
        Q_L[i] = controle_onoff(erro_L, Q_L[i-1] if i > 0 else 0)

        erro_HAG = SP - T1_ModeloHAG[i-1] if i > 0 else erro_real
        Q_HAG[i] = controle_onoff(erro_HAG, Q_HAG[i-1] if i > 0 else 0)

        erro_Smith = SP - T1_ModeloSmith[i-1] if i > 0 else erro_real
        Q_Smith[i] = controle_onoff(erro_Smith, Q_Smith[i-1] if i > 0 else 0)

        erro_Sundaresan = SP - T1_ModeloSundaresan[i-1] if i > 0 else erro_real
        Q_Sundaresan[i] = controle_onoff(erro_Sundaresan, Q_Sundaresan[i-1] if i > 0 else 0)

        erro_ZN = SP - T1_ModeloZN[i-1] if i > 0 else erro_real
        Q_ZN[i] = controle_onoff(erro_ZN, Q_ZN[i-1] if i > 0 else 0)

        # Modelo Não Linear
        if i < L:
            T1_ModeloNL[i] = T1[0]
        else:
            Taquecedor = T1_ModeloNL[i-1] + 273.15
            dTdt = ((Alpha / (m * Cp)) * Q_NL[i]) + ((U * A) / (m * Cp)) * (Tambiente - Taquecedor) \
                   + ((Emissividade * Boltzmann * A) / (m * Cp)) * (Tambiente**4 - Taquecedor**4)
            T1_ModeloNL[i] = T1_ModeloNL[i-1] + dTdt * dt

        # Modelo Linear
        if i <= 0:
            T1_ModeloL[i] = T1[0]
        else:
            T1_ModeloL[i] = T1_ModeloL[i-1] + 0.00499 * np.exp(-tm[i - L] / 123) * Q_L[i] * np.heaviside(i - L, 0)

        # Modelo ZN
        if i <= 0:
            T1_ModeloZN[i] = T1[0]
        else:
            T1_ModeloZN[i] = T1_ModeloZN[i-1] + 0.00272 * np.exp(-tm[i - L] / 217) * Q_ZN[i] * np.heaviside(i - L, 0)

        # Modelo Hägglund
        if i <= 0:
            T1_ModeloHAG[i] = T1[0]
        else:
            T1_ModeloHAG[i] = T1_ModeloHAG[i-1] + 0.00328 * np.exp(-tm[i - L] / 180) * Q_HAG[i] * np.heaviside(i - L, 0)

        # Modelo Smith
        if i <= 0:
            T1_ModeloSmith[i] = T1[0]
        else:
            T1_ModeloSmith[i] = T1_ModeloSmith[i-1] + 0.00340 * np.exp(-tm[i - 25] / 174) * Q_Smith[i] * np.heaviside(i - 25, 0)

        # Modelo Sundaresan
        if i <= 0:
            T1_ModeloSundaresan[i] = T1[0]
        else:
            T1_ModeloSundaresan[i] = T1_ModeloSundaresan[i-1] + 0.00297 * np.exp(-tm[i - 16] / 199) * Q_Sundaresan[i] * np.heaviside(i - 16, 0)

        # Gráficos
        plt.clf()
        ax1 = plt.subplot(2, 1, 1)
        ax1.grid()
        ax1.plot(tm[:i+1], T1[:i+1], 'r-', label='Temperatura Real', alpha=1, zorder=10)
        ax1.plot(tm[:i+1], T1_ModeloNL[:i+1], 'b--', label='Modelo Não Linear', alpha=0.8, zorder=1)
        ax1.plot(tm[:i+1], T1_ModeloL[:i+1], 'g-.', label='Modelo Linear', alpha=0.8, zorder=1)
        ax1.plot(tm[:i+1], T1_ModeloZN[:i+1], 'k-.', label='Modelo ZN', alpha=0.8, zorder=1)
        ax1.plot(tm[:i+1], T1_ModeloHAG[:i+1], 'c-.', label='Modelo HAG', alpha=0.8, zorder=1)
        ax1.plot(tm[:i+1], T1_ModeloSmith[:i+1], 'm-.', label='Modelo Smith', alpha=0.8, zorder=1)
        ax1.plot(tm[:i+1], T1_ModeloSundaresan[:i+1], 'y-.', label='Modelo Sundaresan', alpha=0.8, zorder=1)
        ax1.axhline(SP, color='gray', linestyle='--', label='Setpoint', alpha=1, zorder=2)
        ax1.set_ylabel('Temperatura (°C)')
        leg1 = ax1.legend(loc='best', fontsize=8, ncol=2)
        leg1.set_zorder(100)  # Garante que fique por cima

        ax2 = plt.subplot(2, 1, 2)
        ax2.grid()
        ax2.plot(tm[:i+1], Q1[:i+1], 'r-', label='Potência Aplicada', alpha=1, zorder=10)
        ax2.plot(tm[:i+1], Q_NL[:i+1], 'b--', label='Modelo Não Linear', alpha=0.8, zorder=1)
        ax2.plot(tm[:i+1], Q_L[:i+1], 'g-.', label='Modelo Linear', alpha=0.8, zorder=1)
        ax2.plot(tm[:i+1], Q_ZN[:i+1], 'k-.', label='Modelo ZN', alpha=0.8, zorder=1)
        ax2.plot(tm[:i+1], Q_HAG[:i+1], 'c-.', label='Modelo HAG', alpha=0.8, zorder=1)
        ax2.plot(tm[:i+1], Q_Smith[:i+1], 'm-.', label='Modelo Smith', alpha=0.8, zorder=1)
        ax2.plot(tm[:i+1], Q_Sundaresan[:i+1], 'y-.', label='Modelo Sundaresan', alpha=0.8, zorder=1)
        ax2.set_ylabel('Potência (%)')
        ax2.set_xlabel('Tempo (s)')
        leg2 = ax2.legend(loc='best', fontsize=8, ncol=2)
        leg2.set_zorder(100)  # Garante que fique por cima

        plt.tight_layout()
        plt.draw()
        plt.pause(0.05)

    # Erros médios
    erro_medio_NL = np.mean(np.abs(T1 - T1_ModeloNL))
    erro_medio_L = np.mean(np.abs(T1 - T1_ModeloL))
    erro_medio_HAG = np.mean(np.abs(T1 - T1_ModeloHAG))
    erro_medio_Smith = np.mean(np.abs(T1 - T1_ModeloSmith))
    erro_medio_Sundaresan = np.mean(np.abs(T1 - T1_ModeloSundaresan))
    erro_medio_ZN = np.mean(np.abs(T1 - T1_ModeloZN))
    print(f'\nErro Médio Absoluto - Modelo Não Linear: {erro_medio_NL:.2f} °C')
    print(f'Erro Médio Absoluto - Modelo Linear: {erro_medio_L:.2f} °C')
    print(f'Erro Médio Absoluto - Modelo Hägglund: {erro_medio_HAG:.2f} °C')
    print(f'Erro Médio Absoluto - Modelo Smith: {erro_medio_Smith:.2f} °C')
    print(f'Erro Médio Absoluto - Modelo Sundaresan: {erro_medio_Sundaresan:.2f} °C')
    print(f'Erro Médio Absoluto - Modelo ZN: {erro_medio_ZN:.2f} °C')

    # Salvar gráfico e dados
    plt.savefig('Grafico_ON_OFF_Todos_Modelos.png', dpi=300, bbox_inches='tight')
    np.savetxt('Dados_ON_OFF_Modelos.txt', np.column_stack((tm, T1, T1_ModeloNL, T1_ModeloL, T1_ModeloHAG, T1_ModeloSmith, T1_ModeloSundaresan, T1_ModeloZN, Q1)),
               header='Tempo(s)\tTreal\tTNL\tTL\tTHAG\tTSmith\tTSundaresan\tTZN\tQ1', fmt='%.2f', delimiter='\t')

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
