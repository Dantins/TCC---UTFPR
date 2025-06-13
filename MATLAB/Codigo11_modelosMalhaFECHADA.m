clear; clc;

%% Inicialização do hardware
run('tclab.m');  % Arquivo com as funções T1C(), h1(), led()
% Garante o desligamento seguro do TCLab mesmo em caso de erro ou interrupção (Ctrl+C)
cleanupObj = onCleanup(@() desligarTCLab());

%% Parâmetros de Simulação
run_time = 20;                        % Tempo total da simulação (minutos)
loops = 60 * run_time;                % Número de ciclos (1 ciclo por segundo)
dt = 1.0;                             % Intervalo de tempo em segundos

%% Inicialização dos vetores de dados
tm = zeros(loops,1);                  % Tempo (s)
T1 = zeros(loops,1);                  % Temperatura Real (°C)
T1_ModeloNL = zeros(loops,1);         % Modelo Não Linear (°C)
T1_ModeloL = zeros(loops,1);          % Modelo Linear (°C)
T1_ModeloHAG = zeros(loops,1);        % Modelo Hägglund (°C)
T1_ModeloSmith = zeros(loops,1);      % Modelo Smith (°C)
T1_ModeloSundaresan = zeros(loops,1); % Modelo Sundaresan (°C)
T1_ModeloZN = zeros(loops,1);         % Modelo Ziegler-Nichols (°C)

Q1 = zeros(loops,1);                  % Potência Real Aplicada (%)
Q_NL = zeros(loops,1);                % Potência do Modelo Não Linear (%)
Q_L = zeros(loops,1);                 % Potência do Modelo Linear (%)
Q_HAG = zeros(loops,1);               % Potência do Modelo Hägglund (%)
Q_Smith = zeros(loops,1);             % Potência do Modelo Smith (%)
Q_Sundaresan = zeros(loops,1);        % Potência do Modelo Sundaresan (%)
Q_ZN = zeros(loops,1);                % Potência do Modelo Ziegler-Nichols (%)

%% Setpoint e parâmetros PID
SP = 50;               % Setpoint (°C)
Kp = 1; Ki = 0; Kd = 0; % Parâmetros PID
integral = 0.0;
erro_anterior = 0.0;

%% Parâmetros físicos para o modelo não linear
Alpha = 0.01;
Cp = 500;
A = 0.0012;
m = 0.004;
U = 8;
Emissividade = 0.9;
Boltzmann = 5.67e-8;
Tambiente = 301.15;   % Temperatura ambiente (Kelvin)
L = 15;               % Atraso de transporte (número de ciclos)

%% Criar gráfico interativo
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
        
        %% Leitura de temperatura real
        T1(i) = T1C();
        
        %% Cálculo dos erros
        if i == 1
            Erro_real = SP - T1(i);
            Erro_NL = Erro_real;
            Erro_L = Erro_real;
            Erro_HAG = Erro_real;
            Erro_Smith = Erro_real;
            Erro_Sundaresan = Erro_real;
            Erro_ZN = Erro_real;
        else
            Erro_real = SP - T1(i);
            Erro_NL = SP - T1_ModeloNL(i-1);
            Erro_L = SP - T1_ModeloL(i-1);
            Erro_HAG = SP - T1_ModeloHAG(i-1);
            Erro_Smith = SP - T1_ModeloSmith(i-1);
            Erro_Sundaresan = SP - T1_ModeloSundaresan(i-1);
            Erro_ZN = SP - T1_ModeloZN(i-1);
        end
        
        %% Controle PID para o sistema real
        integral = integral + Erro_real * dt;
        derivada = (Erro_real - erro_anterior) / dt;
        Q_real = Kp * Erro_real + Ki * integral + Kd * derivada;
        erro_anterior = Erro_real;
        Q1(i) = min(max(Q_real, 0), 100);
        h1(Q1(i));  % Aplica a potência no hardware
        
        %% Controle PID para os modelos
        Q_NL(i) = min(max(Kp * Erro_NL, 0), 100);
        Q_L(i) = min(max(Kp * Erro_L, 0), 100);
        Q_HAG(i) = min(max(Kp * Erro_HAG, 0), 100);
        Q_Smith(i) = min(max(Kp * Erro_Smith, 0), 100);
        Q_Sundaresan(i) = min(max(Kp * Erro_Sundaresan, 0), 100);
        Q_ZN(i) = min(max(Kp * Erro_ZN, 0), 100);
        
        %% Modelo Não Linear
        if i < 16
            T1_ModeloNL(i) = T1(1);
        else
            Taquecedor = T1_ModeloNL(i-1) + 273.15; % Converter para Kelvin
            dTdt = ((Alpha/(m*Cp)) * Q_NL(i)) + ((U*A)/(m*Cp)) * (Tambiente - Taquecedor) + ...
                   ((Emissividade*Boltzmann*A)/(m*Cp)) * (Tambiente^4 - Taquecedor^4);
            T1_ModeloNL(i) = T1_ModeloNL(i-1) + dTdt * dt;
        end
        
        %% Modelo Linear
        if i == 1
            T1_ModeloL(i) = T1(1);
        elseif i <= 15
            T1_ModeloL(i) = T1_ModeloL(i-1);
        else
            T1_ModeloL(i) = T1_ModeloL(i-1) + 0.00499 * exp(-tm(i-15)/123) * Q_L(i) * heaviside(i-15 - eps);
        end
        
        %% Modelo Ziegler-Nichols (ZN)
        if i == 1
            T1_ModeloZN(i) = T1(1);
        elseif i <= 15
            T1_ModeloZN(i) = T1_ModeloZN(i-1);
        else
            T1_ModeloZN(i) = T1_ModeloZN(i-1) + 0.00272 * exp(-tm(i-15)/217) * Q_ZN(i) * heaviside(i-15 - eps);
        end
        
        %% Modelo Hägglund
        if i == 1
            T1_ModeloHAG(i) = T1(1);
        elseif i <= 15
            T1_ModeloHAG(i) = T1_ModeloHAG(i-1);
        else
            T1_ModeloHAG(i) = T1_ModeloHAG(i-1) + 0.00328 * exp(-tm(i-15)/180) * Q_HAG(i) * heaviside(i-15 - eps);
        end
        
        %% Modelo Smith
        if i == 1
            T1_ModeloSmith(i) = T1(1);
        elseif i <= 25
            T1_ModeloSmith(i) = T1_ModeloSmith(i-1);
        else
            T1_ModeloSmith(i) = T1_ModeloSmith(i-1) + 0.00340 * exp(-tm(i-25)/174) * Q_Smith(i) * heaviside(i-25 - eps);
        end
        
        %% Modelo Sundaresan
        if i == 1
            T1_ModeloSundaresan(i) = T1(1);
        elseif i <= 16
            T1_ModeloSundaresan(i) = T1_ModeloSundaresan(i-1);
        else
            T1_ModeloSundaresan(i) = T1_ModeloSundaresan(i-1) + 0.00297 * exp(-tm(i-16)/199) * Q_Sundaresan(i) * heaviside(i-16 - eps);
        end
        
        %% Atualização gráfica
        clf;
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
    end
    
    %% Cálculo dos erros médios
    erro_medio_NL = mean(abs(T1 - T1_ModeloNL));
    erro_medio_L = mean(abs(T1 - T1_ModeloL));
    erro_medio_HAG = mean(abs(T1 - T1_ModeloHAG));
    erro_medio_Smith = mean(abs(T1 - T1_ModeloSmith));
    erro_medio_Sundaresan = mean(abs(T1 - T1_ModeloSundaresan));
    erro_medio_ZN = mean(abs(T1 - T1_ModeloZN));
    
    fprintf('\nErro Médio Absoluto entre Temperatura Real e Modelo Não Linear: %.2f °C\n', erro_medio_NL);
    fprintf('Erro Médio Absoluto entre Temperatura Real e Modelo Linear: %.2f °C\n', erro_medio_L);
    fprintf('Erro Médio Absoluto entre Temperatura Real e Modelo Hägglund: %.2f °C\n', erro_medio_HAG);
    fprintf('Erro Médio Absoluto entre Temperatura Real e Modelo Smith: %.2f °C\n', erro_medio_Smith);
    fprintf('Erro Médio Absoluto entre Temperatura Real e Modelo Sundaresan: %.2f °C\n', erro_medio_Sundaresan);
    fprintf('Erro Médio Absoluto entre Temperatura Real e Modelo Ziegler_Nichols: %.2f °C\n', erro_medio_ZN);
    
    % Salvar gráfico final
    saveas(gcf, 'Grafico_TCLab_Todos_Modelos_Comparados_Malha_Fechada.png');
    
    % Salvar dados em arquivo TXT
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
