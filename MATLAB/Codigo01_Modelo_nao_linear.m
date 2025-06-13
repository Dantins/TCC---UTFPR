clear; clc;

%% Inicialização do hardware
run('tclab.m'); % Arquivo com funções T1C(), h1(), led()

% Garantir desligamento ao final, mesmo com Ctrl+C ou erro
cleanupObj = onCleanup(@() desligarTCLab());

%% Parâmetros de Simulação
run_time = 20;                 % Tempo total da simulação (min)
loops = 60 * run_time;         % Número de amostras (1 Hz)
dt = 1.0;                      % Intervalo de tempo (s)

% Inicialização de variáveis
tm = zeros(loops,1);           % Vetor de tempo (s)
T1 = zeros(loops,1);           % Temperatura real (°C)
T1_ModeloNL = zeros(loops,1);  % Temperatura estimada (modelo NL)
Q1 = ones(loops,1) * 50;       % Potência aplicada (%)

% Aplicar potência inicial
h1(Q1(1));

%% Parâmetros físicos do modelo
Tambiente = 301.15;     % Temperatura ambiente (K)
Alpha = 0.01;           % Eficiência do aquecedor
Cp = 500;               % Calor específico (J/kg.K)
A = 0.0012;             % Área de troca térmica (m²)
m = 0.004;              % Massa (kg)
U = 10;                 % Coef. convecção (W/m²K)
Emissividade = 0.9;     % Emissividade
Boltzmann = 5.67e-8;    % Constante de Stefan-Boltzmann

%% Criação do gráfico
figure('Position', [100 100 800 600]);

%% Loop principal de aquisição e simulação
start_time = tic;
prev_time = toc(start_time);

for i = 1:loops
    % Espera até completar 1 segundo
    sleep_time = max(0.01, 1.0 - (toc(start_time) - prev_time));
    pause(sleep_time);

    % Atualizar tempo
    t = toc(start_time);
    prev_time = t;
    tm(i) = t;

    % Leitura da temperatura real
    T1(i) = T1C();

    % Atualização do modelo não linear
    if i == 1
        T1_ModeloNL(i) = T1(i);
    else
        Taquecedor = T1_ModeloNL(i-1) + 273.15; % Converter para Kelvin
        Q = Q1(i);
        
        % Equação diferencial do modelo
        dTdt = ((Alpha / (m * Cp)) * Q) ...
             + ((U * A) / (m * Cp)) * (Tambiente - Taquecedor) ...
             + ((Emissividade * Boltzmann * A) / (m * Cp)) * (Tambiente^4 - Taquecedor^4);
        
        % Integração de Euler
        T1_ModeloNL(i) = T1_ModeloNL(i-1) + dTdt * dt;
    end

    % Atualização gráfica
    clf;
    subplot(2,1,1);
    plot(tm(1:i), T1(1:i), 'r-', 'LineWidth', 2); hold on;
    plot(tm(1:i), T1_ModeloNL(1:i), 'b--', 'LineWidth', 2);
    ylabel('Temperatura (°C)', 'FontSize', 14);
    legend('Temperatura Real', 'Modelo Não Linear', 'FontSize', 12);
    grid on;

    subplot(2,1,2);
    plot(tm(1:i), Q1(1:i), 'r-', 'LineWidth', 2);
    ylabel('Potência (%)', 'FontSize', 14);
    xlabel('Tempo (s)', 'FontSize', 14);
    legend('Potência (%)', 'FontSize', 12);
    grid on;

    drawnow;
end

%% Pós-processamento
erro_medio = mean(abs(T1 - T1_ModeloNL));
fprintf('\nErro Médio Absoluto entre Temperatura Real e Modelo: %.2f °C\n', erro_medio);

% Salvar gráfico
saveas(gcf, 'Grafico_TCLab_ModeloNL.png');

% Salvar dados em TXT
dados = [tm, T1, T1_ModeloNL, Q1];
header = {'Tempo_s', 'Temperatura_Real_C', 'Temperatura_Modelo_C', 'Potencia_pct'};
T = array2table(dados, 'VariableNames', header);
writetable(T, 'Dados_simulacao_ModeloNL.txt', 'Delimiter', '\t');

disp('Simulação concluída.');

%% Função de desligamento seguro
function desligarTCLab()
    try
        h1(0);
        led(0);
        disp('Dispositivos TCLab desligados com segurança.');
    catch
        disp('Erro ao tentar desligar os dispositivos.');
    end
end
