clear; clc;

%% Inicialização do hardware
run('tclab.m');  % Arquivo com as funções T1C(), h1(), led()
% Garante o desligamento seguro do TCLab mesmo em caso de erro ou interrupção (Ctrl+C)
cleanupObj = onCleanup(@() desligarTCLab());

%% Parâmetros de Simulação
run_time = 20;                        % Tempo total da simulação (minutos)
loops = 60 * run_time;                % Número de ciclos (1 ciclo por segundo)
dt = 1.0;                             % Intervalo de tempo em segundos

% Vetores de armazenamento
tm = zeros(loops,1);                  % Tempo (s)
T1 = zeros(loops,1);                  % Temperatura real (°C)
T1_ModeloNL = zeros(loops,1);         % Modelo Não Linear (°C)
T1_ModeloL = zeros(loops,1);          % Modelo Linear (°C)
T1_ModeloHAG = zeros(loops,1);        % Modelo Hägglund (°C)
T1_ModeloSmith = zeros(loops,1);      % Modelo Smith (°C)
T1_ModeloSundaresan = zeros(loops,1); % Modelo Sundaresan (°C)
T1_ModeloZN = zeros(loops,1);         % Modelo Ziegler-Nichols (°C)
Q1 = ones(loops,1) * 50;              % Potência fixa: 50%

% Aplicar potência inicial ao hardware e aos modelos
h1(Q1(1));
Q = 50;

% Parâmetros do modelo não linear
Tambiente = 301.15;  % Temperatura ambiente (29 °C em Kelvin)
Alpha = 0.01;
Cp = 500;
A = 0.0012;
m = 0.004;
U = 8;
Emissividade = 0.9;
Boltzmann = 5.67e-8;

%% Criação do Gráfico
figure('Position',[100 100 1200 800]);
drawnow;

%% Loop principal de aquisição e modelagem
start_time = tic;
prev_time = toc(start_time);

try
    for i = 1:loops
        % Aguarda para manter o intervalo dt
        sleep_time = max(0.01, dt - (toc(start_time) - prev_time));
        pause(sleep_time);
        
        % Atualizar tempo
        t = toc(start_time);
        prev_time = t;
        tm(i) = t;
        
        % Ler temperatura real do TCLab
        T1(i) = T1C();
        
        %% Modelo Não Linear
        if i < 16  % Para i = 1 a 15
            T1_ModeloNL(i) = T1(1);
        else
            Taquecedor = T1_ModeloNL(i-1) + 273.15;  % Converter para Kelvin
            dTdt = ((Alpha/(m*Cp)) * Q) + ((U*A)/(m*Cp)) * (Tambiente - Taquecedor) + ...
                   ((Emissividade*Boltzmann*A)/(m*Cp)) * (Tambiente^4 - Taquecedor^4);
            T1_ModeloNL(i) = T1_ModeloNL(i-1) + dTdt * dt;
        end
        
        %% Modelo Linear
        if i == 1
            T1_ModeloL(i) = T1(1);  % Condição inicial
        elseif i <= 15
            T1_ModeloL(i) = T1_ModeloL(i-1);
        else
            T1_ModeloL(i) = T1_ModeloL(i-1) + 0.00499 * exp(-tm(i-15)/123) * Q * heaviside(i-15 - eps);
        end
        
        %% Modelo Ziegler-Nichols (ZN)
        if i == 1
            T1_ModeloZN(i) = T1(1);  % Condição inicial
        elseif i <= 15
            T1_ModeloZN(i) = T1_ModeloZN(i-1);
        else
            T1_ModeloZN(i) = T1_ModeloZN(i-1) + 0.00272 * exp(-tm(i-15)/217) * Q * heaviside(i-15 - eps);
        end
        
        %% Modelo Hägglund
        if i == 1
            T1_ModeloHAG(i) = T1(1);  % Condição inicial
        elseif i <= 15
            T1_ModeloHAG(i) = T1_ModeloHAG(i-1);
        else
            T1_ModeloHAG(i) = T1_ModeloHAG(i-1) + 0.00328 * exp(-tm(i-15)/180) * Q * heaviside(i-15 - eps);
        end
        
        %% Modelo Smith
        if i == 1
            T1_ModeloSmith(i) = T1(1);  % Condição inicial
        elseif i <= 25
            T1_ModeloSmith(i) = T1_ModeloSmith(i-1);
        else
            T1_ModeloSmith(i) = T1_ModeloSmith(i-1) + 0.00340 * exp(-tm(i-25)/174) * Q * heaviside(i-25 - eps);
        end
        
        %% Modelo Sundaresan
        if i == 1
            T1_ModeloSundaresan(i) = T1(1);  % Condição inicial
        elseif i <= 16
            T1_ModeloSundaresan(i) = T1_ModeloSundaresan(i-1);
        else
            T1_ModeloSundaresan(i) = T1_ModeloSundaresan(i-1) + 0.00297 * exp(-tm(i-16)/199) * Q * heaviside(i-16 - eps);
        end
        
        %% Atualização gráfica
        clf;
        subplot(2,1,1);
        plot(tm(1:i), T1(1:i), 'r-', 'LineWidth', 2); hold on;
        plot(tm(1:i), T1_ModeloNL(1:i), 'b--', 'LineWidth', 2);
        plot(tm(1:i), T1_ModeloL(1:i), 'g--', 'LineWidth', 2);
        plot(tm(1:i), T1_ModeloSmith(1:i), 'y--', 'LineWidth', 2);
        plot(tm(1:i), T1_ModeloHAG(1:i), 'c--', 'LineWidth', 2);
        plot(tm(1:i), T1_ModeloZN(1:i), 'm--', 'LineWidth', 2);
        plot(tm(1:i), T1_ModeloSundaresan(1:i), 'k--', 'LineWidth', 2);
        ylabel('Temperatura (°C)', 'FontSize', 14);
        legend('Temperatura Real', 'Modelo Não Linear', 'Modelo Linear', 'Modelo Smith', ...
               'Modelo Hägglund', 'Modelo ZN', 'Modelo Sundaresan', 'Location', 'best', 'FontSize', 12);
        grid on;
        
        subplot(2,1,2);
        plot(tm(1:i), Q1(1:i), 'r-', 'LineWidth', 2);
        ylabel('Potência (%)', 'FontSize', 14);
        xlabel('Tempo (s)', 'FontSize', 14);
        legend('Potência (%)', 'Location', 'best', 'FontSize', 12);
        grid on;
        
        drawnow;
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
    fprintf('Erro Médio Absoluto entre Temperatura Real e Modelo Ziegler_Nichols: %.2f °C\n', erro_medio_ZN);
    fprintf('Erro Médio Absoluto entre Temperatura Real e Modelo Hägglund: %.2f °C\n', erro_medio_HAG);
    fprintf('Erro Médio Absoluto entre Temperatura Real e Modelo Smith: %.2f °C\n', erro_medio_Smith);
    fprintf('Erro Médio Absoluto entre Temperatura Real e Modelo Sundaresan: %.2f °C\n', erro_medio_Sundaresan);
    
    % Salvar gráfico em PNG
    saveas(gcf, 'Grafico_TCLab_Todos_Modelos_Comparados_Malha_Aberta.png');
    
    % Salvar dados em arquivo TXT
    dados = [tm, T1, T1_ModeloNL, T1_ModeloL, T1_ModeloHAG, T1_ModeloSmith, T1_ModeloSundaresan, T1_ModeloZN, Q1];
    header = {'Tempo_s', 'Temperatura_Real_C', 'Modelo_Nao_Linear_C', 'Modelo_Linear_C', 'Modelo_Hagglund_C', ...
              'Modelo_Smith_C', 'Modelo_Sundaresan_C', 'Modelo_ZN_C', 'Potencia_pct'};
    T_table = array2table(dados, 'VariableNames', header);
    writetable(T_table, 'Dados_simulacao_10_Modelos_Comparados_Malha_Aberta.txt', 'Delimiter', '\t');
    
catch ME
    disp(['Erro: ' ME.message]);
    rethrow(ME);
end

disp('Simulação concluída.');

%% Função de desligamento seguro do TCLab
function desligarTCLab()
    try
        h1(0);    % Desliga o aquecedor
        led(0);   % Apaga o LED
        disp('Dispositivos TCLab desligados com segurança.');
    catch
        disp('Erro ao tentar desligar os dispositivos.');
    end
end
