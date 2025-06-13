clear; clc;

%% Inicialização do hardware
run('tclab.m');  % Arquivo com as funções T1C(), h1(), led()
% Garante o desligamento seguro do TCLab mesmo em caso de erro ou interrupção (Ctrl+C)
cleanupObj = onCleanup(@() desligarTCLab());

%% Parâmetros de Simulação
run_time = 20;                % Tempo total da simulação (minutos)
loops = 60 * run_time;        % Número de ciclos (1 ciclo por segundo)
dt = 1;                       % Intervalo de tempo (segundos)

%% Inicialização dos vetores de dados
tm                      = zeros(loops,1);  % Tempo
T1                      = zeros(loops,1);  % Temperatura Real
T1_ModeloNL             = zeros(loops,1);  % Modelo Não Linear
T1_ModeloL              = zeros(loops,1);  % Modelo Linear
T1_ModeloHAG            = zeros(loops,1);  % Modelo Hägglund
T1_ModeloSmith          = zeros(loops,1);  % Modelo Smith
T1_ModeloSundaresan     = zeros(loops,1);  % Modelo Sundaresan
T1_ModeloZN             = zeros(loops,1);  % Modelo Ziegler-Nichols

Q1                      = zeros(loops,1);  % Potência Real Aplicada
Q_NL                    = zeros(loops,1);  % Potência para Modelo Não Linear
Q_L                     = zeros(loops,1);  % Potência para Modelo Linear
Q_HAG                   = zeros(loops,1);  % Potência para Modelo Hägglund
Q_Smith                 = zeros(loops,1);  % Potência para Modelo Smith
Q_Sundaresan            = zeros(loops,1);  % Potência para Modelo Sundaresan
Q_ZN                    = zeros(loops,1);  % Potência para Modelo ZN

%% Setpoint e parâmetros PID ideal
SP = 50;        % Setpoint (°C)
Kp = 20;        % Ganho Proporcional
Ti = 50;        % Tempo Integral
Td = 10;        % Tempo Derivativo
integral = 0.0;
erro_anterior = 0.0;

% Variáveis PID para os modelos
integral_NL         = 0.0;
integral_L          = 0.0;
integral_HAG        = 0.0;
integral_Smith      = 0.0;
integral_Sundaresan = 0.0;
integral_ZN         = 0.0;
erro_anterior_NL         = 0.0;
erro_anterior_L          = 0.0;
erro_anterior_HAG        = 0.0;
erro_anterior_Smith      = 0.0;
erro_anterior_Sundaresan = 0.0;
erro_anterior_ZN         = 0.0;

%% Parâmetros físicos para o modelo não linear
Alpha = 0.01;
Cp = 500;
A  = 0.0012;
m  = 0.004;
U  = 8;
Emissividade = 0.9;
Boltzmann = 5.67e-8;
Tambiente = 301.15;   % Kelvin
L = 15;               % Atraso de transporte (em número de ciclos)

%% Criação do gráfico interativo
figure('Position',[100 100 1000 700]);
drawnow;

%% Início da simulação
start_time = tic;
prev_time = toc(start_time);

try
    for i = 1:loops
        %% Controle de tempo
        sleep_time = max(0.01, 1.0 - (toc(start_time) - prev_time));
        pause(sleep_time);
        t = toc(start_time);
        prev_time = t;
        tm(i) = t;
        
        %% Leitura da temperatura real do dispositivo
        T1(i) = T1C();
        
        %% Cálculo dos erros (PID)
        if i == 1
            Erro_real       = SP - T1(i);
            Erro_NL         = Erro_real;
            Erro_L          = Erro_real;
            Erro_HAG        = Erro_real;
            Erro_Smith      = Erro_real;
            Erro_Sundaresan = Erro_real;
            Erro_ZN         = Erro_real;
        else
            Erro_real       = SP - T1(i);
            Erro_NL         = SP - T1_ModeloNL(i-1);
            Erro_L          = SP - T1_ModeloL(i-1);
            Erro_HAG        = SP - T1_ModeloHAG(i-1);
            Erro_Smith      = SP - T1_ModeloSmith(i-1);
            Erro_Sundaresan = SP - T1_ModeloSundaresan(i-1);
            Erro_ZN         = SP - T1_ModeloZN(i-1);
        end
        
        %% Controle PID ideal para o sistema real
        integral = integral + Erro_real * dt;
        derivada = (Erro_real - erro_anterior) / dt;
        Q_real = Kp * (Erro_real + (1/Ti)*integral + Td * derivada);
        % Ajuste para evitar saturação
        if Q_real > 100 && Erro_real > 0
            integral = integral - Erro_real * dt;
        elseif Q_real < 0 && Erro_real < 0
            integral = integral - Erro_real * dt;
        end
        Q1(i) = min(max(Q_real, 0), 100);
        erro_anterior = Erro_real;
        h1(Q1(i));  % Aplica a potência real
        
        %% Controle PID ideal para os modelos (utilizando a função local pid_model_ideal)
        [Q_NL(i), integral_NL, erro_anterior_NL] = pid_model_ideal(Erro_NL, erro_anterior_NL, integral_NL, Kp, Ti, Td);
        [Q_L(i),  integral_L,  erro_anterior_L]  = pid_model_ideal(Erro_L,  erro_anterior_L,  integral_L,  Kp, Ti, Td);
        [Q_HAG(i), integral_HAG, erro_anterior_HAG] = pid_model_ideal(Erro_HAG, erro_anterior_HAG, integral_HAG, Kp, Ti, Td);
        [Q_Smith(i), integral_Smith, erro_anterior_Smith] = pid_model_ideal(Erro_Smith, erro_anterior_Smith, integral_Smith, Kp, Ti, Td);
        [Q_Sundaresan(i), integral_Sundaresan, erro_anterior_Sundaresan] = pid_model_ideal(Erro_Sundaresan, erro_anterior_Sundaresan, integral_Sundaresan, Kp, Ti, Td);
        [Q_ZN(i), integral_ZN, erro_anterior_ZN] = pid_model_ideal(Erro_ZN, erro_anterior_ZN, integral_ZN, Kp, Ti, Td);
        
        %% Atualização dos modelos
        % Modelo Não Linear
        if i < L
            T1_ModeloNL(i) = T1(1);
        else
            Taquecedor = T1_ModeloNL(i-1) + 273.15;
            dTdt = ((Alpha/(m*Cp)) * Q_NL(i)) + ((U*A)/(m*Cp)) * (Tambiente - Taquecedor) + ...
                   ((Emissividade*Boltzmann*A)/(m*Cp)) * (Tambiente^4 - Taquecedor^4);
            T1_ModeloNL(i) = T1_ModeloNL(i-1) + dTdt * dt;
        end
        
        % Modelos Linear, ZN, Hägglund, Smith e Sundaresan
        if i == 1
            T1_ModeloL(i) = T1(i);
            T1_ModeloZN(i) = T1(i);
            T1_ModeloHAG(i) = T1(i);
            T1_ModeloSmith(i) = T1(i);
            T1_ModeloSundaresan(i) = T1(i);
        else
            if i <= L
                T1_ModeloL(i) = T1_ModeloL(i-1);
                T1_ModeloZN(i) = T1_ModeloZN(i-1);
                T1_ModeloHAG(i) = T1_ModeloHAG(i-1);
            else
                T1_ModeloL(i) = T1_ModeloL(i-1) + 0.00499 * exp(-tm(i-L)/123) * Q_L(i) * heaviside(i - L - eps);
                T1_ModeloZN(i) = T1_ModeloZN(i-1) + 0.00272 * exp(-tm(i-L)/217) * Q_ZN(i) * heaviside(i - L - eps);
                T1_ModeloHAG(i) = T1_ModeloHAG(i-1) + 0.00328 * exp(-tm(i-L)/180) * Q_HAG(i) * heaviside(i - L - eps);
            end
            if i <= 25
                T1_ModeloSmith(i) = T1_ModeloSmith(i-1);
            else
                T1_ModeloSmith(i) = T1_ModeloSmith(i-1) + 0.00340 * exp(-tm(i-25)/174) * Q_Smith(i) * heaviside(i - 25 - eps);
            end
            if i <= 16
                T1_ModeloSundaresan(i) = T1_ModeloSundaresan(i-1);
            else
                T1_ModeloSundaresan(i) = T1_ModeloSundaresan(i-1) + 0.00297 * exp(-tm(i-16)/199) * Q_Sundaresan(i) * heaviside(i - 16 - eps);
            end
        end
        
        %% Atualização dos gráficos
        clf;
        % Gráfico da temperatura
        subplot(2,1,1);
        plot(tm(1:i), T1(1:i), 'r-', 'LineWidth', 2); hold on;
        plot(tm(1:i), T1_ModeloNL(1:i), 'b--', 'LineWidth', 2);
        plot(tm(1:i), T1_ModeloL(1:i), 'g-.', 'LineWidth', 2);
        plot(tm(1:i), T1_ModeloSmith(1:i), 'y-.', 'LineWidth', 2);
        plot(tm(1:i), T1_ModeloHAG(1:i), 'c-.', 'LineWidth', 2);
        plot(tm(1:i), T1_ModeloZN(1:i), 'm-.', 'LineWidth', 2);
        plot(tm(1:i), T1_ModeloSundaresan(1:i), 'k-.', 'LineWidth', 2);
        yline(SP, 'Color', '#808080', 'LineStyle', '--', 'LineWidth', 1.5, 'DisplayName', 'Setpoint');
        ylabel('Temperatura (°C)', 'FontSize', 14);
        legend('Temperatura Real', 'Modelo Não Linear', 'Modelo Linear', 'Modelo Smith', ...
               'Modelo Hägglund', 'Modelo ZN', 'Modelo Sundaresan', 'Location', 'best', 'FontSize', 12);
        grid on;
        
        % Gráfico da potência
        subplot(2,1,2);
        plot(tm(1:i), Q1(1:i), 'r-', 'LineWidth', 2); hold on;
        plot(tm(1:i), Q_NL(1:i), 'b--', 'LineWidth', 2);
        plot(tm(1:i), Q_L(1:i), 'g-.', 'LineWidth', 2);
        plot(tm(1:i), Q_ZN(1:i), 'm-.', 'LineWidth', 2);
        plot(tm(1:i), Q_HAG(1:i), 'c-.', 'LineWidth', 2);
        plot(tm(1:i), Q_Smith(1:i), 'y-.', 'LineWidth', 2);
        plot(tm(1:i), Q_Sundaresan(1:i), 'k-.', 'LineWidth', 2);
        ylabel('Potência (%)', 'FontSize', 14);
        xlabel('Tempo (s)', 'FontSize', 14);
        legend('Potência Real', 'Modelo Não Linear', 'Modelo Linear', 'Modelo ZN', 'Modelo Hägglund', ...
               'Modelo Smith', 'Modelo Sundaresan', 'Location', 'best', 'FontSize', 12);
        grid on;
        
        drawnow;
        pause(0.05);
    end  % fim do laço for

    %% Cálculo dos erros médios
    erro_medio_NL         = mean(abs(T1 - T1_ModeloNL));
    erro_medio_L          = mean(abs(T1 - T1_ModeloL));
    erro_medio_HAG        = mean(abs(T1 - T1_ModeloHAG));
    erro_medio_Smith      = mean(abs(T1 - T1_ModeloSmith));
    erro_medio_Sundaresan = mean(abs(T1 - T1_ModeloSundaresan));
    erro_medio_ZN         = mean(abs(T1 - T1_ModeloZN));
    
    fprintf('\nErro Médio Absoluto entre Temperatura Real e Modelo Não Linear: %.2f °C\n', erro_medio_NL);
    fprintf('Erro Médio Absoluto entre Temperatura Real e Modelo Linear: %.2f °C\n', erro_medio_L);
    fprintf('Erro Médio Absoluto entre Temperatura Real e Modelo Hägglund: %.2f °C\n', erro_medio_HAG);
    fprintf('Erro Médio Absoluto entre Temperatura Real e Modelo Smith: %.2f °C\n', erro_medio_Smith);
    fprintf('Erro Médio Absoluto entre Temperatura Real e Modelo Sundaresan: %.2f °C\n', erro_medio_Sundaresan);
    fprintf('Erro Médio Absoluto entre Temperatura Real e Modelo Ziegler_Nichols: %.2f °C\n', erro_medio_ZN);
    
    % Salvar o gráfico final
    saveas(gcf, 'Grafico_TCLab_Todos_Modelos_Comparados_Malha_Fechada_PID.png');
    
    % Salvar os dados em arquivo TXT
    dados = [tm, T1, T1_ModeloNL, T1_ModeloL, T1_ModeloHAG, T1_ModeloSmith, T1_ModeloSundaresan, T1_ModeloZN, Q1];
    header = {'Tempo_s', 'Temperatura_Real_C', 'Modelo_Nao_Linear_C', 'Modelo_Linear_C', 'Modelo_Hagglund_C', ...
              'Modelo_Smith_C', 'Modelo_Sundaresan_C', 'Modelo_ZN_C', 'Potencia_pct'};
    T_table = array2table(dados, 'VariableNames', header);
    writetable(T_table, 'Dados_simulacao_11_Modelos_Comparados_Malha_Fechada.txt', 'Delimiter', '\t');
    
catch ME
    disp(['Erro: ' ME.message]);
    rethrow(ME);
end

disp('Simulação concluída.');

%% Função de desligamento seguro do TCLab
function desligarTCLab()
    try
        h1(0);
        led(0);
        disp('Dispositivos TCLab desligados com segurança.');
    catch
        disp('Erro ao tentar desligar os dispositivos.');
    end
end

%% Função PID ideal para modelos
function [Q, integral_acum, erro_novo] = pid_model_ideal(Erro, erro_ant, integral_acum, Kp, Ti, Td)
    dt = 1;
    integral_acum = integral_acum + Erro * dt;
    derivada = (Erro - erro_ant) / dt;
    Q = Kp * (Erro + (1/Ti)*integral_acum + Td*derivada);
    if Q > 100 && Erro > 0
        integral_acum = integral_acum - Erro*dt;
    elseif Q < 0 && Erro < 0
        integral_acum = integral_acum - Erro*dt;
    end
    Q = min(max(Q,0),100);
    erro_novo = Erro;
end
