import tclab
import numpy as np
import time
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes, mark_inset

# ------------------- CONTROLE COM TCLAB -------------------
with tclab.TCLab() as a:
    run_time = 20  # Tempo em minutos
    loops = int(60.0 * run_time)

    tm = np.zeros(loops)
    T1 = np.zeros(loops)
    Q1 = np.zeros(loops)

    P_term = np.zeros(loops)
    I_term = np.zeros(loops)
    D_term = np.zeros(loops)

    SP = 50
    Kp = 17.36
    Ti = 30
    Td = 7.5

    integral = 0.0
    erro_anterior = 0.0
    dt = 1

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

            erro = SP - T1[i]

            P = Kp * erro
            integral += erro * dt
            I = Kp * (1 / Ti) * integral
            derivada = (erro - erro_anterior) / dt
            D = Kp * Td * derivada

            P_term[i] = P
            I_term[i] = I
            D_term[i] = D

            Q = P + I + D

            # Anti-windup
            if Q > 100 and erro > 0:
                integral -= erro * dt
            elif Q < 0 and erro < 0:
                integral -= erro * dt

            Q = np.clip(Q, 0, 100)
            Q1[i] = Q

            erro_anterior = erro
            a.Q1(Q)

            # Plot em tempo real
            plt.clf()
            ax1 = plt.subplot(2, 1, 1)
            ax1.grid()
            ax1.plot(tm[:i+1], T1[:i+1], 'r-', label='Temperatura Real')
            ax1.axhline(SP, color='gray', linestyle='--', label='Setpoint')
            ax1.set_ylabel('Temperatura (°C)')
            ax1.legend(loc='best')

            ax2 = plt.subplot(2, 1, 2)
            ax2.grid()
            ax2.plot(tm[:i+1], Q1[:i+1], 'r-', label='Potência Aplicada')
            ax2.set_ylabel('Potência (%)')
            ax2.set_xlabel('Tempo (s)')
            ax2.legend(loc='best')

            plt.draw()
            plt.pause(0.05)

        plt.ioff()

    except KeyboardInterrupt:
        print('Interrompido pelo usuário.')
    except Exception as e:
        print(f'Erro: {e}')
        raise

print('Simulação encerrada. Hardware desligado com segurança.')

# ---------------- GRÁFICO COM ZOOM AJUSTADO ----------------
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 7))

# Gráfico da temperatura
ax1.grid()
ax1.plot(tm, T1, 'r-', label='Temperatura Real')
ax1.axhline(SP, color='gray', linestyle='--', label='Setpoint')
ax1.set_ylabel('Temperatura (°C)')
ax1.legend(loc='best')

# Zoom no início, agora reposicionado mais acima
axins = inset_axes(ax1, width="40%", height="30%",
                   loc='upper right',
                   bbox_to_anchor=(0.6, 0.4, 0.4, 0.4),
                   bbox_transform=ax1.transAxes, borderpad=2)

axins.plot(tm, T1, 'r-')
axins.axhline(SP, color='gray', linestyle='--')
axins.set_xlim(0, 30)
axins.set_ylim(min(T1[tm <= 30]) - 1, max(T1[tm <= 30]) + 1)
axins.grid()
mark_inset(ax1, axins, loc1=3, loc2=4, fc="none", ec="0.5")

# Gráfico da potência
ax2.grid()
ax2.plot(tm, Q1, 'r-', label='Potência Aplicada')
ax2.set_ylabel('Potência (%)')
ax2.set_xlabel('Tempo (s)')
ax2.legend(loc='best')

plt.tight_layout()
plt.savefig('Grafico_TCLab_PID_ZN_COM_ZOOM.png', dpi=300)
plt.show()

# ---------------- GRÁFICO DOS TERMOS PID (DIVIDIDO) ----------------
fig, axs = plt.subplots(4, 1, figsize=(10, 10), sharex=True)

# Proporcional
axs[0].plot(tm, P_term, 'g-', label='Proporcional (P)')
axs[0].set_ylabel('Ação (%)')
axs[0].set_ylim(-50, 400)
axs[0].grid()
axs[0].legend(loc='best')

# Integral
axs[1].plot(tm, I_term, 'b-', label='Integral (I)')
axs[1].set_ylabel('Ação (%)')
axs[1].set_ylim(-50, 150)
axs[1].grid()
axs[1].legend(loc='best')

# Derivativo
axs[2].plot(tm, D_term, color='olive', label='Derivativo (D)')
axs[2].set_ylabel('Ação (%)')
axs[2].set_ylim(-50, 150)
axs[2].grid()
axs[2].legend(loc='best')

# PID Total
axs[3].plot(tm, Q1, 'r-', label='PID (Saída)')
axs[3].set_ylabel('Ação (%)')
axs[3].set_xlabel('Tempo (s)')
axs[3].set_ylim(-50, 150)
axs[3].grid()
axs[3].legend(loc='best')

fig.suptitle('Componentes do Controle PID', fontsize=14)
plt.tight_layout(rect=[0, 0, 1, 0.96])
plt.savefig('Grafico_Componentes_PID_DIVIDIDOS.png', dpi=300)
plt.show()

# ---------------- SALVANDO DADOS ----------------
np.savetxt('Dados_simulacao_PID_ZN_ZOOM.txt',
           np.column_stack((tm, T1, Q1, P_term, I_term, D_term)),
           header='Tempo(s)\tTemperatura(°C)\tPotência(%)\tP\tI\tD',
           fmt='%.2f', delimiter='\t')
