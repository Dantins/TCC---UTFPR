clear; clc;

%% Inicialização do hardware
run('tclab.m'); % Arquivo com as funções T1C(), h1(), led()
% Garante o desligamento seguro do TCLab mesmo em caso de erro ou interrupção (Ctrl+C)
cleanupObj = onCleanup(@() desligarTCLab());

%% Parâmetros de Simulação
run_time = 20;                 % Tempo total da simulação (min)
loops = 60 * run_time;         % Número de amostras (1 Hz)
dt = 1.0;                      % Intervalo de tempo (s)

% Inicialização de variáveis
tm = zeros(loops,1);           % Vetor de tempo (s)
T1 = zeros(loops,1);           % Temperatura real (°C)
T1_ModeloNL = zeros(loops,1);  % Temperatura modelada não linear (°C)
T1_ModeloL = zeros(loops,1);   % Temperatura modelada linear (°C)
Q1 = ones(loops,1) * 50;       % Potência fixa aplicada (%)
u = zeros(loops,1);            % Entrada no modelo linear (se necessário)

% Aplicar potência inicial ao hardware e aos modelos
h1(Q1(1));
Q = 50;

%% Parâmetros do Modelo Não Linear
Tambiente = 301.15;  % Temperatura ambiente (K) – 25 °C em Kelvin
Alpha = 0.01;
Cp = 500;
A = 0.0012;
m = 0.004;
U = 8;
Emissividade = 0.9;
Boltzmann = 5.67e-8;

%% Criação do Gráfico
figure('Position',[100 100 800 600]);

%% Loop principal de aquisição e simulação
start_time = tic;
prev_time = toc(start_time);

for i = 1:loops
    % Tempo de espera até completar dt segundos
    sleep_time = max(0.01, dt - (toc(start_time) - prev_time));
    pause(sleep_time);
    
    % Atualizar tempo
    t = toc(start_time);
    prev_time = t;
    tm(i) = t;
    
    % Leitura da temperatura real do TCLab
    T1(i) = T1C();
    
    % --- Modelo Não Linear ---
    % Para os primeiros 16 ciclos (equivalente a índices 0 a 15 em Python) utiliza condição inicial
    if i <= 16  
        T1_ModeloNL(i) = T1(1);
    else
        Taquecedor = T1_ModeloNL(i-1) + 273.15;  % Converter para Kelvin
        dTdt = ((Alpha / (m * Cp)) * Q) + ((U * A) / (m * Cp)) * (Tambiente - Taquecedor) + ...
                ((Emissividade * Boltzmann * A) / (m * Cp)) * (Tambiente^4 - Taquecedor^4);
        T1_ModeloNL(i) = T1_ModeloNL(i-1) + dTdt * dt;
    end
    
    % --- Modelo Linear ---
    % Condição inicial no primeiro instante
    if i == 1
        T1_ModeloL(i) = T1(1);
    elseif i <= 15
        % Enquanto o atraso não ocorre, mantém o valor anterior
        T1_ModeloL(i) = T1_ModeloL(i-1);
    else
        % Atualização somente após 15 ciclos; a função heaviside do MATLAB
        T1_ModeloL(i) = T1_ModeloL(i-1) + 0.00499 * exp(-tm(i-15)/123) * Q * heaviside(i-15);
    end
    
    % --- Atualização Gráfica ---
    clf;
    subplot(2,1,1);
    plot(tm(1:i), T1(1:i), 'r-', 'LineWidth', 2); hold on;
    plot(tm(1:i), T1_ModeloNL(1:i), 'b--', 'LineWidth', 2);
    plot(tm(1:i), T1_ModeloL(1:i), 'g--', 'LineWidth', 2);
    ylabel('Temperatura (°C)', 'FontSize', 14);
    legend('Temperatura Real', 'Modelo Não Linear', 'Modelo Linear', 'FontSize', 12);
    grid on;
    
    subplot(2,1,2);
    plot(tm(1:i), Q1(1:i), 'r-', 'LineWidth', 2);
    xlabel('Tempo (s)', 'FontSize', 14);
    ylabel('Potência (%)', 'FontSize', 14);
    legend('Potência (%)', 'FontSize', 12);
    grid on;
    
    drawnow;
end

%% Pós-processamento
erro_medio_NL = mean(abs(T1 - T1_ModeloNL));
erro_medio_L  = mean(abs(T1 - T1_ModeloL));
fprintf('\nErro Médio Absoluto entre Temperatura Real e Modelo Não Linear Com Atraso: %.2f °C\n', erro_medio_NL);
fprintf('Erro Médio Absoluto entre Temperatura Real e Modelo Linear Com Atraso: %.2f °C\n', erro_medio_L);

% Salvar gráfico em PNG
saveas(gcf, 'Grafico_TCLab_Modelo_Linear_Com_Atraso.png');

% Salvar dados em arquivo TXT
dados = [tm, T1, T1_ModeloNL, T1_ModeloL, Q1];
header = {'Tempo_s', 'Temperatura_Real_C', 'Modelo_Nao_Linear_C', 'Modelo_Linear_C', 'Potencia_pct'};
T = array2table(dados, 'VariableNames', header);
writetable(T, 'Dados_simulacao_04.txt', 'Delimiter', '\t');

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
