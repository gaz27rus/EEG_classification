D1 = table2array(D1); %Перевод таблицы данных ЭЭГ в массив Matlab
D2 = table2array(D2); % (имя переменной поменять в зависимости от имени импротируемого файла)
data = {D1, D2}; %Транспонирование матрицы данных ЭЭГ к виду 
% (число_каналов х количество отсчётов), если требуется

EEG_Data = prepare_data(data);
%% 
time_parameters = struct('start_bias', 5, ... %отступ от начала файла в с 
                         'trial_bias', [0, 1, 2], ... %смещение для нарезки проб в с
                         'trial_time', 2, ... %длительность пробы в с
                         'labels', [3,6], ...
                         'srate', 250); %Частота дискретизации

EEG_Metki = prepare_metki(EEG_Data, time_parameters);%создаются массивы отсчётов для всех смещений
%% 
%Определение параметров
%set(0,'DefaultFigureVisible','on');%Отключить отображение графиков в матлабе
resolution = '.png'; %Разрешение файла
single_save = 0;% 1 - Раздельное сохранение каналов, 0 - слитное сохранение

EEG_parameters = struct('start_time', 0, ... %Начало интересующего промежутка (в секундах) от начала пробы
                        'end_time', 2.4, ... %Конец интересующего промежутка (в секундах) от начала пробы
                        'srate', 250, ... %Частота дискретизации
                        'n_trials', length(EEG_Metki), ... %Количество проб
                        'n_channels', 19); %Число каналов
EEG_parameters.n = EEG_parameters.end_time*EEG_parameters.srate - EEG_parameters.start_time*EEG_parameters.srate; ... %Число отсчётов из требуемого промежутка
%%
%Блок нарезки исходных данных на требуемые пробы)
output_data = cell(1, length(EEG_Data));
for i = 1 : length(output_data)
    output_data{i} = cell(1, length(time_parameters.trial_bias));
end

for i = 1 : length(EEG_Metki)
    temp_sig = cell2mat(EEG_Data(i));
    t = EEG_Metki{i};
    X = cell(1, length(time_parameters.trial_bias));
    for j = 1 : length(t)
        metki = t{j};
        X{j} = crop_data(temp_sig, metki, EEG_parameters);
    end
    output_data{i} = X;
end

%На выходе получаем массив
%(число_каналов х число_проб х число_отсчётов_в_пробе)
%% 
final_data = cell(1, length(output_data));
data_size = 30; %сколько данных нужно взять из каждого набора

for i = 1 : length(output_data)
    temp = output_data{i};
    temp_data = zeros(EEG_parameters.n_channels, data_size*length(output_data{i}), EEG_parameters.n);
    for j = 1 : length(temp)
        x = temp{j};
        temp_data(:, ((j-1)*data_size+1):(j*data_size), :) = x(:, 1:data_size, :);
    end
    final_data{i} = temp_data;
end

%% Вайвлеты
channels = [5, 15]; %Задаём массив каналов, из которых требуется извлечь ЧВ преобразование
global merged_img;

for k = 1 : length(final_data)
    temp = final_data{k};

    for i = 1 : length(temp(1, :, 1))
        temp_path = 'C:\Users\Глеб\Desktop\Wavelet\Data\images\temp\'; %Директория для сохранения промежуточных изображений (потом их можно удалить)
        if ~exist(temp_path, 'dir')
            mkdir(temp_path);
        end
        
        resize_size = [302, 412];
        merged_size = [resize_size(1)*length(channels), resize_size(2)];
            
        merged_img = zeros(merged_size(1), merged_size(2), 3, 'uint8');
    
        for j = 1 : length(channels)
            image = tiledlayout(1, 1);
            sig = squeeze(temp(j,i,:));
        
            t = 0:1/EEG_parameters.srate:(length(sig)-1)*1/EEG_parameters.srate; %временной отрезок отрисовки
            [cfs,f] = cwt(sig,EEG_parameters.srate,'amor', VoicesPerOctave = 40);
                
            t_s = 0.5; %начало временного промежутка отрисовки ЧВ-графика
            t_e = 2.2; %конец временного промежутка отрисовки ЧВ-графика
            f_s = 3; %нижняя частотота отрисовки
            f_e = 30; %верхняя частотота отрисовки
        
            t_range = t(t >= t_s & t <= t_e);
            f_range = f(f >= f_s & f <= f_e);
        
            cfs_range = cfs(f >= f_s & f <= f_e, t >= t_s & t <= t_e);
        
            ax = nexttile();
            ax.Position = [60 10 216 300];
            imagesc(t_range, f_range, abs(cfs_range));
            caxis([0 8]); %диапазон амплитуд
            set(gca,'ydir','normal', 'YScale','log');
            xlabel('Time (s)');
            ylabel('Frequency (Hz)');
            colorbar;
            truesize([216 300]);
        
        
            temp_file_name = sprintf('%s%s%s%s%s%s%s', num2str(k), 'set_',  ...
                'file', num2str(i), '_', num2str(j), resolution);
            out_file=[temp_path, temp_file_name];
            exportgraphics(image, out_file);
        
            crop_and_save(time_parameters.labels(k),i,j,temp_path,temp_file_name, channels, ...
                single_save);
        end
    end
end
%% Сохранение файлов для матлаба
channels_output = [5, 15];
test = cell(1, length(EEG_Data));

for k = 1 : length(final_data)
    X = final_data{k};
    out = zeros(length(X(1, :, 1)), length(channels_output)*length(X(1, 1, :)));

    for i = 1 : length(X(1, :, 1))
        for j = 1 : length(channels_output)
            temp = squeeze(X(channels_output(j), i, :)).';
            out(i, length(temp)*(j-1)+1:length(temp)*j) = temp;
        end
    end
    
    test{k} = out;

    temp_path = 'C:\Users\Глеб\Desktop\Wavelet\Data\images\txt_data_for_classifier\'; %Директория для сохранения промежуточных изображений (потом их можно удалить)
    if ~exist(temp_path, 'dir')
        mkdir(temp_path);
    end
    file_name = sprintf('%s%s%s%s%s%s%s', temp_path, num2str(time_parameters.labels(k)), ...
        '_', num2str(length(X(1, :, 1))), 'trial_', num2str(length(temp)), ...
        'samples');
    writematrix(out, sprintf('%s%s',file_name, '.txt'), 'Delimiter','tab'); % write the matrix to a text file
    trial_samples = EEG_parameters.n;
    save(sprintf('%s%s',file_name, '.mat'), 'X', ...
        'EEG_parameters');
end
%% 
%% Сохранение файлов в общем формате
channels_output = [5, 15];
test1 = cell(1, length(EEG_Data));

for k = 1 : length(final_data)
    X = final_data{k};
    out = zeros(length(X(1, :, 1))*length(channels_output), length(X(1, 1, :)));

    for i = 1 : length(channels_output)
        temp = squeeze(X(channels_output(j), :, :));
        out(length(X(1, :, 1))*(i-1)+1:length(X(1, :, 1))*(i), :) = temp;
    end
    
    test1{k} = out;

    temp_path = 'C:\Users\Глеб\Desktop\Wavelet\Data\images\txt_data_for_art\'; %Директория для сохранения промежуточных изображений (потом их можно удалить)
    if ~exist(temp_path, 'dir')
        mkdir(temp_path);
    end
    file_name = sprintf('%s%s%s%s%s%s%s', temp_path, num2str(time_parameters.labels(k)), ...
        '_', num2str(length(X(1, :, 1))), 'trial_', num2str(length(temp)), ...
        'samples');
    writematrix(out, sprintf('%s%s',file_name, '.txt'), 'Delimiter','tab'); % write the matrix to a text file
    trial_samples = EEG_parameters.n;
    save(sprintf('%s%s',file_name, '.mat'), 'X', ...
        'EEG_parameters');
end
%% 
function EEG_Data = prepare_data(input_arrays)

    EEG_Data = cell(1, length(input_arrays));

    for i = 1: length(EEG_Data)
        EEG_Data{i} = input_arrays{i}.';
    end
end
%% 
function EEG_Metki = prepare_metki(input_arrays, time_parameters)

    EEG_Metki = cell(1, length(input_arrays));

    for i = 1 : length(EEG_Metki)
        EEG_Metki{i} = cell(1, length(time_parameters.trial_bias));
    end

    for k = 1 : length(EEG_Metki)
        x = input_arrays{k};
        
        EEG_Metki_trial = cell(1, length(time_parameters.trial_bias));
        for i = 1 : length(time_parameters.trial_bias)
            t = zeros(1, (int32((length(x(1, :))-(time_parameters.start_bias+time_parameters.trial_bias(i))*time_parameters.srate)/(time_parameters.trial_time*time_parameters.srate))-1));
            for j = 1 : length(t)
                t(j) = time_parameters.start_bias+time_parameters.trial_bias(i)+(j-1)*time_parameters.trial_time;
            end
            EEG_Metki_trial{i} = t;
        end
        EEG_Metki{k} = EEG_Metki_trial;
    end
end
%% 
function data = crop_data(in_data, in_metki, EEG_Parameters)

    data = zeros(EEG_Parameters.n_channels,length(in_metki),EEG_Parameters.n);
    for i = 1 : EEG_Parameters.n_channels
        X_trials = zeros(length(in_metki), EEG_Parameters.n);
        for j = 1 : length(in_metki)-1
            var = in_data(i, (int32(in_metki(1,j)*EEG_Parameters.srate+EEG_Parameters.start_time*EEG_Parameters.srate)-1):(int32(in_metki(1,j+1)*EEG_Parameters.srate+EEG_Parameters.end_time*EEG_Parameters.srate)+1));
            X_trials(j,:) = var(1:EEG_Parameters.n);

            %var = in_data((int32(in_metki(1,j)*EEG_Parameters.srate+EEG_Parameters.start_time*EEG_Parameters.srate)-1):(int32(in_metki(1,j)*EEG_Parameters.srate+EEG_Parameters.end_time*EEG_Parameters.srate)+1));
            %X_trials(j,:) = var(1:EEG_Parameters.n);
        end
        data(i,:,:) = X_trials;
    end
end
%% 
function crop_and_save(k,i,j,temp_path,temp_file_name, channels, ...
    single_save)
    img = imread(sprintf('%s%s', temp_path, temp_file_name));
  
    global merged_img;
    resize_size = [301, 411];
    cropped_img = imcrop(img, [68, 10, 410, 300]);%Обрезание картинки

    if single_save == 1%Если сохранение раздельное  
        file_name = sprintf('%s%s%s%s%s%s', num2str(k), ...
            '_', num2str(i), 'tr_', num2str(channels(j)), 'ch.png');%Имя файла
        
        path = sprintf('%s%s%s', 'C:\Users\Глеб\Desktop\Wavelet\Data\images\', num2str(k), '\');%Директория сохранения с делением на классы
        
        if ~exist(path, 'dir')%Проверка, существует ли директория
            mkdir(path);
        end
        
        imwrite(cropped_img, sprintf('%s%s', path, file_name));%Сохранение файла
    else%Если слитное сохранение
        canvas_x = (j-1) * resize_size(1) + 1;
        merged_img(canvas_x:canvas_x+resize_size(1)-1, 1:resize_size(2), :) = cropped_img;

        if j == length(channels)%Если последний канал - объединяем изображение и сохраняем
            file_name = sprintf('%s%s%s%s', num2str(k), '_', ...
                num2str(i), 'tr.png');

            path = sprintf('%s%s%s', 'C:\Users\Глеб\Desktop\Wavelet\Data\images\', num2str(k), '\');%Директория сохранения с делением на классы

            if ~exist(path, 'dir')
                mkdir(path);
            end

            imwrite(merged_img, sprintf('%s%s', path, file_name));
        end
    end 
end