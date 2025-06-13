clear; clc;

%% Inicialização do hardware
run('tclab.m');  % Arquivo com as funções T1C(), h1(), led()
% Garante o desligamento seguro do TCLab mesmo em caso de erro ou interrupção (Ctrl+C)
cleanupObj = onCleanup(@() desligarTCLab());

%% Parâmetros de Simulação
run_time = 20;                 % Tempo total da simulação (min)
loops = 60 * run_time;         % Número de ciclos (1 ciclo por segundo)
dt = 1.0;                      % Intervalo de tempo (s)

% Inicialização das variáveis
tm = zeros(loops,1);           % Vetor para armazenar o tempo (s)
T1 = zeros(loops,1);           % Temperatura real (°C)
T1_ModeloHAG = zeros(loops,1);   % Temperatura do modelo Hägglund (°C)
Q1 = ones(loops,1) * 50;       % Potência fixa em 50%
% u permanece não utilizado neste exemplo

% Aplicar potência inicial ao hardware
h1(Q1(1));
Q = 50;  % Potência aplicada aos modelos

%% Criação do Gráfico
figure('Position', [100 100 800 600]);
drawnow;

%% Loop principal de aquisição e simulação
start_time = tic;
prev_time = toc(start_time);

try
    for i = 1:loops
        % Tempo de espera para manter o intervalo de dt segundos
        sleep_time = max(0.01, dt - (toc(start_time) - prev_time));
        pause(sleep_time);
        
        % Atualizar tempo
        t = toc(start_time);
        prev_time = t;
        tm(i) = t;
        
        % Ler temperatura real do TCLab
        T1(i) = T1C();
        
        % --- Modelo Hägglund ---
        if i == 1
            % Condição inicial
            T1_ModeloHAG(i) = T1(1);
        elseif i <= 15
            % Enquanto o atraso não ocorre, mantém o valor anterior
            T1_ModeloHAG(i) = T1_ModeloHAG(i-1);
        else
            % Para i > 15: atualiza o modelo utilizando o fator de atraso.
            % A função heaviside é ajustada com "-eps" para que heaviside(0-eps) retorne 0.
            T1_ModeloHAG(i) = T1_ModeloHAG(i-1) + 0.00328 * exp(-tm(i-15)/180) * Q * heaviside(i-15 - eps);
        end
        
        % --- Atualização Gráfica ---
        clf;
        subplot(2,1,1);
        plot(tm(1:i), T1(1:i), 'r-', 'LineWidth', 2); hold on;
        plot(tm(1:i), T1_ModeloHAG(1:i), 'g--', 'LineWidth', 2);
        ylabel('Temperatura (°C)', 'FontSize', 14);
        legend('Temperatura Real', 'Modelo Hägglund', 'FontSize', 12);
        grid on;
        
        subplot(2,1,2);
        plot(tm(1:i), Q1(1:i), 'r-', 'LineWidth', 2);
        ylabel('Potência (%)', 'FontSize', 14);
        xlabel('Tempo (s)', 'FontSize', 14);
        legend('Potência (%)', 'FontSize', 12);
        grid on;
        
        drawnow;
    end
    
    % --- Cálculo do Erro e Salvamento ---
    erro_medio_HAG = mean(abs(T1 - T1_ModeloHAG));
    fprintf('Erro Médio Absoluto entre Temperatura Real e Modelo Hägglund: %.2f °C\n', erro_medio_HAG);
    
    % Salvar gráfico em PNG
    saveas(gcf, 'Grafico_TCLab_Modelo_HAG.png');
    
    % Salvar dados em arquivo TXT
    dados = [tm, T1, T1_ModeloHAG, Q1];
    header = {'Tempo_s', 'Temperatura_Real_C', 'Modelo_HAG_C', 'Potencia_pct'};
    T_data = array2table(dados, 'VariableNames', header);
    writetable(T_data, 'Dados_simulacao_07_HAG.txt', 'Delimiter', '\t');
    
catch ME
    disp(['Erro: ', ME.message]);
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
