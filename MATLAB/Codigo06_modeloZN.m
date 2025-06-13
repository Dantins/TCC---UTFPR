clear; clc;

%% Inicialização do hardware
run('tclab.m');  % Arquivo com as funções T1C(), h1(), led()
% Garante o desligamento seguro do TCLab mesmo em caso de erro ou interrupção (Ctrl+C)
cleanupObj = onCleanup(@() desligarTCLab());

%% Parâmetros de Simulação
run_time = 20;                 % Tempo total da simulação (min)
loops = 60 * run_time;         % Número de ciclos (1 ciclo por segundo)
dt = 1.0;                      % Intervalo de tempo em segundos

% Inicialização das variáveis
tm = zeros(loops,1);           % Vetor para armazenar o tempo (s)
T1 = zeros(loops,1);           % Temperatura real (°C)
T1_ModeloZN = zeros(loops,1);   % Temperatura do modelo Ziegler-Nichols (°C)
Q1 = ones(loops,1) * 50;       % Potência fixa em 50%
% (u permanece não utilizado neste exemplo)
  
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
        % Tempo de espera para que o intervalo seja de dt segundos
        sleep_time = max(0.01, dt - (toc(start_time) - prev_time));
        pause(sleep_time);
        
        % Atualizar tempo
        t = toc(start_time);
        prev_time = t;
        tm(i) = t;
        
        % Ler temperatura real do TCLab
        T1(i) = T1C();
        
        % --- Modelo Ziegler-Nichols (ZN) ---
        if i == 1
            % Condição inicial
            T1_ModeloZN(i) = T1(1);
        elseif i < 16
            % Para i < 16 (ou seja, enquanto o atraso não ocorre), manter valor anterior
            T1_ModeloZN(i) = T1_ModeloZN(i-1);
        else
            % A atualização ocorre somente após 15 ciclos; como i-15 > 0, heaviside(i-15) retorna 1
            T1_ModeloZN(i) = T1_ModeloZN(i-1) + 0.00272 * exp(-tm(i-15) / 217) * Q * heaviside(i-15);
        end
        
        % --- Atualizar Gráfico ---
        clf;
        subplot(2,1,1);
        plot(tm(1:i), T1(1:i), 'r-', 'LineWidth', 2); hold on;
        plot(tm(1:i), T1_ModeloZN(1:i), 'g--', 'LineWidth', 2);
        ylabel('Temperatura (°C)', 'FontSize', 14);
        legend('Temperatura Real', 'Modelo Ziegler\_Nichols', 'FontSize', 12);
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
    erro_medio_ZN = mean(abs(T1 - T1_ModeloZN));
    fprintf('Erro Médio Absoluto entre Temperatura Real e Modelo Ziegler_Nichols: %.2f °C\n', erro_medio_ZN);
    
    % Salvar gráfico em PNG
    saveas(gcf, 'Grafico_TCLab_Modelo_ZN.png');
    
    % Salvar dados em arquivo TXT
    dados = [tm, T1, T1_ModeloZN, Q1];
    header = {'Tempo_s', 'Temperatura_Real_C', 'Modelo_ZN_C', 'Potencia_pct'};
    T_data = array2table(dados, 'VariableNames', header);
    writetable(T_data, 'Dados_simulacao_06_ZN.txt', 'Delimiter', '\t');
    
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
