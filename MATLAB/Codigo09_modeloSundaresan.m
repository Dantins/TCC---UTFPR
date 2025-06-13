clear; clc;

%% Inicialização do hardware
run('tclab.m');  % Arquivo com as funções T1C(), h1(), led()
% Garante o desligamento seguro do TCLab mesmo em caso de erro ou interrupção (Ctrl+C)
cleanupObj = onCleanup(@() desligarTCLab());

%% Parâmetros de Simulação
run_time = 20;                        % Tempo total da simulação (min)
loops = 60 * run_time;                % Número de ciclos (1 ciclo por segundo)
dt = 1.0;                             % Intervalo de tempo em segundos

% Inicialização das variáveis
tm = zeros(loops,1);                  % Vetor para armazenar o tempo (s)
T1 = zeros(loops,1);                  % Temperatura real (°C)
T1_ModeloSundaresan = zeros(loops,1); % Temperatura do modelo Sundaresan (°C)
Q1 = ones(loops,1) * 50;              % Potência fixa em 50%
% (u não é utilizado neste exemplo)

% Aplicar potência inicial ao hardware
h1(Q1(1));
Q = 50;  % Potência aplicada aos modelos

% Parâmetro do modelo Sundaresan (não utilizado diretamente no cálculo,
% mas pode servir para documentar o atraso em número de ciclos, ex.: 18)
Atraso_Sundaresan = 18;

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
        
        % --- Modelo Sundaresan ---
        % Em Python (indexando de 0): se i == 0, usa a condição inicial;
        % Para 0 < i < 16 => np.heaviside(i-16,0) retorna 0; para i >= 16, o modelo é atualizado.
        % Em MATLAB (indexando de 1): se i == 1, condição inicial; para i <= 16, manter o valor anterior;
        % para i >= 17, atualizar o modelo.
        if i == 1
            T1_ModeloSundaresan(i) = T1(1);  % condição inicial
        elseif i <= 16
            T1_ModeloSundaresan(i) = T1_ModeloSundaresan(i-1);
        else
            % Note: tm(i-16) corresponde a tm[i-16] em Python, considerando o deslocamento de índice.
            T1_ModeloSundaresan(i) = T1_ModeloSundaresan(i-1) + ...
                0.00297 * exp(-tm(i-16) / 199) * Q;
        end
        
        % --- Atualização Gráfica ---
        clf;
        subplot(2,1,1);
        plot(tm(1:i), T1(1:i), 'r-', 'LineWidth', 2); hold on;
        plot(tm(1:i), T1_ModeloSundaresan(1:i), 'g--', 'LineWidth', 2);
        ylabel('Temperatura (°C)', 'FontSize', 14);
        legend('Temperatura Real', 'Modelo Sundaresan', 'FontSize', 12);
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
    erro_medio_Sundaresan = mean(abs(T1 - T1_ModeloSundaresan));
    fprintf('Erro Médio Absoluto entre Temperatura Real e Modelo Sundaresan: %.2f °C\n', erro_medio_Sundaresan);
    
    % Salvar gráfico em PNG
    saveas(gcf, 'Grafico_TCLab_Modelo_Sundaresan.png');
    
    % Salvar dados em arquivo TXT
    dados = [tm, T1, T1_ModeloSundaresan, Q1];
    header = {'Tempo_s', 'Temperatura_Real_C', 'Modelo_Sundaresan_C', 'Potencia_pct'};
    T_data = array2table(dados, 'VariableNames', header);
    writetable(T_data, 'Dados_simulacao_09_Sundaresan.txt', 'Delimiter', '\t');
    
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
