% 
clear;
% result02FVA = [0.00, 0.83; 1.37, 2.09; 2.60 ,3.57; 4.00, 4.76; 5.33, 6.18; 6.68, 7.17];
% result01MDA = [0.00 ,0.45; 0.81, 1.53; 1.85, 2.69; 2.86, 3.78; 4.15, 4.84; 5.14 ,5.58];
% result03MAB = [0.00	,1.03; 1.42, 2.46; 2.80	,4.21 ;4.52, 6.81 ;7.14,	8.22; 8.50, 9.37];
% result06FTP = [0.00	,1.52; 1.92, 3.91; 4.35,	6.18; 6.60 ,8.67 ;9.14,	10.94; 11.33,12.75];
silenceMDV = [0.00, 0.88; 1.34, 2.35; 2.82, 3.76; 4.13, 5.04; 5.50, 6.41; 6.79, 7.42];
silenceMTT = [0.00, 0.93; 1.42, 2.59; 3.00, 4.71; 5.11, 6.26; 6.66, 8.04; 8.39, 9.27];
silenceFQT = [0.00, 0.46; 0.99, 1.56; 2.13, 2.51; 2.93, 3.79; 4.38, 4.77; 5.22, 5.79];
silenceFTN = [0.00, 0.59; 0.97, 1.76; 2.11, 3.44; 3.77, 4.70; 5.13, 5.96; 6.28, 6.78];

for k = 0: 0
%         voiceDetected('G:\Document5\XLTH\Cuoi-Ki-Thuc-Hanh\Tin-Hieu-Huan-Luyen\01MDA.wav',result01MDA);
%         voiceDetected('G:\Document5\XLTH\Cuoi-Ki-Thuc-Hanh\Tin-Hieu-Huan-Luyen\02FVA.wav',result02FVA);
%         voiceDetected('G:\Document5\XLTH\Cuoi-Ki-Thuc-Hanh\Tin-Hieu-Huan-Luyen\03MAB.wav',result03MAB);
%         voiceDetected('G:\Document5\XLTH\Cuoi-Ki-Thuc-Hanh\Tin-Hieu-Huan-Luyen\06FTB.wav',result06FTP);
    voiceDetected('G:\THKT\30FTN.wav',silenceFTN);
    voiceDetected('G:\THKT\42FQT.wav',silenceFQT);
    voiceDetected('G:\THKT\44MTT.wav',silenceMTT);
    voiceDetected('G:\THKT\45MDV.wav',silenceMDV);
end
%Ham chinh dung de phan biet khoang lang va tieng noi
function voiceDetected(fileName , standardSignal)
[audioIn, Fs] = audioread(fileName);
fileName
% 1.---Cài đặt các giá trị cơ bản---
audioIn = audioIn./abs(max(audioIn));
samples = length(audioIn);
frame_duration = 0.02;
frame_length = round(Fs * frame_duration); % so mau trong 1 frame
frameTotalWithoutFrameShift = floor(samples / frame_length);
frame_total = 2*frameTotalWithoutFrameShift - 1;% tong so frame duoc chia ra
Weight = 1;
%----Thực hiện gọi hàm ----

% 2.---Tính STE của mỗi frame---
STE_PowFrame_Matrix = computeSTE(audioIn, frame_total, frame_length);

% 3.---Chuẩn hóa STE về với biên độ [0, 1]---
STE_PowFrame_Matrix = Standard_function(STE_PowFrame_Matrix);

% 4.---Làm mịn giá trị đặc trưng STE
baseSTE = Compute_BaseSTE(STE_PowFrame_Matrix,frame_length);

% 5.---Threshold---
threshHold =Compute_Threshold(STE_PowFrame_Matrix, Weight);

%6.---Phân tích tiếng nói và khoảng lặng
checkSpeechArray = AnalysisVoice_Function(frame_total , STE_PowFrame_Matrix , threshHold);
% xac dinh cac vi tri (bat dau, ket thuc) cua khoang lang
% Khai báo biến thời gian . Đưa trục Ox về đơn vị thời gian
time = [1 : length(audioIn)] / Fs; % Đưa về đơn vị thời gian ( s )
time_STE = (1 : length(baseSTE))/ Fs; % Đưa về đơn vị thời gian ( s )

% 7. Dựng hình vẽ đồ thị
silenceIndexArray = findSilenceIndex(checkSpeechArray, frame_total);
voiceArray = DetectVoiceIndexArray(silenceIndexArray);
% 8. Gọi hàm tính F0 
[F0_voice ,F0,P1,freqV,freqU,P_unvoice]= compute_F0(voiceArray , fileName);
F0_mean = mean(F0_voice);
F0_std =  std(F0_voice);
figure('name',fileName);

subplot(4,2,[1 , 2]);        % Vẽ hai đồ thị chung một figure
plot(time,audioIn);hold on;
plot(time_STE,baseSTE , 'r','LineWidth',1.5);
title('Short-Time-Energy (STE)');  % Gán nhãn tiêu đề đồ thị
legend('Voice Signal','STE','location','northeast'); % Chú thích các thành phần trên đồ thị
xlabel('Thời gian (s)'); % Gán nhãn đơn vị trục thời gian
ylabel('Biên độ');% Gán nhãn đơn vị trục biên độ
subplot(4,2,[3 , 4]);     % Vẽ hai đồ thị chung một figure
p1=plot(audioIn);
hold on;
%Vẽ biên chuẩn file lab
for j = 1: length(standardSignal)
    p2=xline(standardSignal(j, 1) * Fs, 'r','LineWidth',2);
    xline(standardSignal(j, 2) * Fs, 'r','LineWidth',2);
end
% Vẽ biên theo thuật toán
for j = 1: size(silenceIndexArray)
    start =  silenceIndexArray(j, 1);
    endIndex = silenceIndexArray(j, 2);
    p3=xline(((start - 1) / 2) * frame_length, 'g','LineWidth',2);
    xline(((endIndex * frame_length) - (frame_length * (endIndex - 1) / 2)), 'g','LineWidth',2);
end
hold off;
legend([p1,p2,p3],{'Voice Signal' ,'By Teacher', 'By Student'});  % Chú thích các thành phần trên đồ thị
xlabel('Frame'); % Gán nhãn đơn vị trục thời gian
ylabel('Biên độ'); % Gán nhãn đơn vị trục thời gian
title('Phân đoạn speech và silence (STE)'); % Gán nhãn tiêu đề đồ thị
subplot(4,2,5);
plot(freqV(1:length(freqV)/10), P1(1:length(P1)/10));
findpeaks(P1(1:length(P1)/10), freqV(1:length(freqV)/10), 'NPeaks', 3, 'MinPeakDistance', 110, 'MinPeakHeight', 1);
legend('Spectrum','Peak');
xlabel('Tần số (Hz)'); 
ylabel('Biên độ'); 
title('Phổ biên độ khung Voice');
subplot(4,2,6)
plot(freqU(1:length(freqU)/10), P_unvoice(1:length(P_unvoice)/10));
legend('Spectrum');
xlabel('Tần số (Hz)'); 
ylabel('Biên độ'); 
title('Phổ biên độ khung Unvoice');
subplot(4,2,[7,8]);
stem(F0,'.','LineStyle', 'none', 'MarkerFace', 'b');
hold on;
xlabel('Frame'); 
ylabel('Tần số (Hz)'); 
legend('F0');
title(['Tìm F0 trên miền tần số (Hz), F0mean = ',num2str(F0_mean),' F0std = ',num2str(F0_std)]);
% figure(2);
% [histSTE, x_STE] = hist(STE_PowFrame_Matrix, round(length(STE_PowFrame_Matrix)/0.15)); % Tần suất xuất hiện ( hist STE ) giá trị STE mỗi frame
% bar(histSTE,200);
% ylabel('Hist Value'); 
% legend('Hist');
% title('Hist Graph');
end

%----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%                                                                                                             HÀM XỬ LÝ PHÂN TÍCH                                                                                                            %                                           
%----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------%
% 1. Hàm xác định khoảng lặng
function silenceIndexArray = findSilenceIndex(checkSpeechArray, frameTotal)
silenceIndexArray = [];
indexSilence = 1;
stepForSkip = 0;
for i = 1 : frameTotal
    if(stepForSkip > 0)
        stepForSkip = stepForSkip - 1;
        continue;
    end
    % Kiem tra xem khoang lang >= 300ms hay khong
    if(checkSpeechArray(i) == 0)
        count = i;
        while(count < frameTotal && checkSpeechArray(count + 1) == 0)
            count = count + 1;
        end
        if(count - i >= 14)
            silenceIndexArray(indexSilence, 1) = i;
            silenceIndexArray(indexSilence, 2) = count;
            indexSilence = indexSilence + 1;
            stepForSkip = count - i;
        end
    end
end
end
% 2. Hàm tính STE  mỗi frame
function STE_PowFrame_Matrix = computeSTE(x, frameTotal, frameLength)
STE_PowFrame_Matrix = zeros(1, frameTotal); % tinh nang luong cua moi frame
for i = 1 : frameTotal
    startIndex = (frameLength * (i - 1) / 2) + 1;
    endIndex =  startIndex + frameLength - 1 ;
    frameI = x(startIndex : endIndex);
    % tien hanh tinh STE:
    STE_PowFrame_Matrix(i) = sum(frameI.^2);
end
end
% 3. Hàm tính ngưỡng (Threshold) sử dụng thuật toán Histogram
function [threshHold] = Compute_Threshold(STE_PowFrame_Matrix, Weight)
[histSTE, x_STE] = hist(STE_PowFrame_Matrix, round(length(STE_PowFrame_Matrix)/1.5)); % Tần suất xuất hiện ( hist STE ) giá trị STE mỗi frame

% tại các vị trí x_STE.
% vecto histSTE : lưu tần suất xuất hiện ( số lần xuất hiện ) giá trị STE.
% của mỗi frame ( STE_PowFrame_Matrix) tại vị trí x_STE ( vecto ).
maximaHistSTE1 = 0;
maximaHistSTE2 = 0;
maximaIndex1 = 0; % Vị trí cực đại cục bộ thứ 1
maximaIndex2 = 0; % Vị trí cực đại cục bộ thứ 2
%Tìm cực đại cục bộ thứ nhất và thứ hai nằm cùng 1 frame
for i = 2 : length(histSTE) - 1  % Duyệt kết quả đồ thị tần suất ( histSTE)
    previous = i - 1;
    next = i + 1;
    while(histSTE(i) == histSTE(next)) % Xét vị trí histSTE thứ i và histSTE liền kề
        next = next + 1;
    end
    if(histSTE(i) > histSTE(previous) && histSTE(i) > histSTE(next)) % Kiểm tra giá trị tại histSTE thứ i so với giá trị tại histSTE trước và sau
        if(maximaIndex1 == 0)
            maximaHistSTE1 = histSTE(i);
            maximaIndex1 = i;
        else
            maximaHistSTE2 = histSTE(i);
            maximaIndex2 = i;
            break;
        end
    end
    i = next;
end
maximaHistSTE1 = x_STE(maximaIndex1); % Kết quả giá trị cực đại cục bộ thứ nhất
maximaHistSTE2 = x_STE(maximaIndex2); % Kết quả giá trị cực đại cục bộ thứ hai
% B2: Áp dụng công thức : T = (W * M1 + M2) / (W + 1)
threshHold = (Weight * maximaHistSTE1 + maximaHistSTE2) / (Weight + 1);
end
% 4. Hàm tính baseSTE 
function [baseSTE] = Compute_BaseSTE(STE_PowFrame_Matrix,frame_length)
    for i = 1: length(STE_PowFrame_Matrix)
        startIndex = (frame_length * (i - 1) / 2) + 1;
        endIndex =  startIndex + frame_length - 1 ;
        baseSTE(startIndex : endIndex) = STE_PowFrame_Matrix(i);
    end
end
function [std] = std1(F0voice)
    sum = 0;
    sd = 0;
    for i = 1 : length(F0voice)
        sum = sum + F0voice(i);
    end
    mean = sum / length(F0voice);
    for i = 1 : length(F0voice)
        sd = sd + (F0voice(i) - mean)^2;
    end
    std = sqrt(sd / length(F0voice));
    std = std-5;
end
% 5. Hàm xác định biên 
function [Point] = Point(Speech, frame_number, frame_duration)  % Mảng các điểm tiếng nói và khoảng lặng
    Point = [];
    j = 1;
    for i = 2:frame_number                                % Duyệt tất cả các khung
        if(Speech(i) == 1 && Speech(i-1)==0)
            Point(2*j - 1)= (i-1)*frame_duration;         % Lưu vị trí biên đầu của tiếng nói
            j = j + 1;
        end
    end
    j = 1;
    for i = 2:frame_number
        if (Speech(i) == 0 && Speech(i-1) == 1)
            Point(2*j) = (i-1)*frame_duration;             % Lưu vị trí biên cuối của tiếng nói
            j = j + 1;
        end
    end
end
% 6. Hàm chuẩn hóa STE về biên [0:1]
function [STE_PowFrame_Matrix] = Standard_function(STE_PowFrame_Matrix)
    minEnergy = min(STE_PowFrame_Matrix); % Giá trị Energy min
    maxEnergy = max(STE_PowFrame_Matrix);% Giá trị Energy max
    for i = 1 : length(STE_PowFrame_Matrix) % Duyệt từng phần tử mảng STE_PowFrame_Matrix
        STE_PowFrame_Matrix(i) = (STE_PowFrame_Matrix(i) - minEnergy) / (maxEnergy - minEnergy);
    end
end
% 7. Hàm đưa ra quyết định tiếng nói hay khoảng lặng => Đưa ra quyết định
% voice hay unvoice
function [checkSpeechArray] =  AnalysisVoice_Function(frame_total , STE_PowFrame_Matrix , threshHold)
    checkSpeechArray = zeros(1, frame_total); % Mảng kiểm tra tiếng nói và khoảng lặng
    for i = 1 : frame_total % Duyệt từng frame
        if(STE_PowFrame_Matrix(i) >threshHold )  % So sánh với giá trị ngưỡng
            checkSpeechArray(i) = 1; % Giá trị STE lớn hơn giá trị ngưỡng => Lưu bằng 1
        else
            checkSpeechArray(i) = 0; % % Giá trị STE lớn hơn giá trị ngưỡng => Lưu bằng 0
        end
    end
end
% 8. Lấy mảng khoảng tiếng nói sử dụng lại kết quả câu a
function [voiceArrayIndexArray] = DetectVoiceIndexArray(silenceIndexArray)
    i = 1;
    for k = 1 : size(silenceIndexArray)-1 
        voiceArrayIndexArray(i,1) = silenceIndexArray(i,2) + 1;
        voiceArrayIndexArray(i,2) = silenceIndexArray(i+1,1)-1;
        i = i + 1;
    end
end
% 9. Kiểm tra vùng là tiếng nói 
function isCheck  = CheckVoice(voiceArray , k)
     isCheck = 0;
     for  i = 1 : size(voiceArray) 
         if k <= voiceArray(i , 2) && k >=  voiceArray(i , 1)
             isCheck = 1;
             break;
         end
     end
end
% 10. Hàm cửa sổ hamming 
function ham = my_haming(frame_length)
    ham = zeros(1,frame_length);
    for i=1:frame_length
        ham(i) = 0.54 - 0.46*cos(2*pi*(i-1)/(frame_length));
    end
end
% 11. Hàm tìm F0 sử dụng FFT 
function [F0_Voice,F0,P_Voice,freqV,freqU,P_Unvoice] = compute_F0(voiceIndexArray , fileName)
    warning('off');
    [audioIn , Fs] = audioread(fileName);
    frame_length = round(0.02*Fs);    %Độ dài khung tín hiệu (20ms)
    frameTotalWithoutFrameShift = floor(length(audioIn) / frame_length);%Tổng số frame được chia ra 
    frame_total = 2*frameTotalWithoutFrameShift - 1;%Xử lý overlap
    hamm = my_haming(frame_length); %Hamming windows
    N = 32768;
    F0_Voice = [];
    F0=[];
    P_Unvoice = [];
    checkSilence = 1;
    index=1;      %Chỉ số của số phần tử trong F0_Voice
    for k = 1 : frame_total  % Duyệt qua tổng frame
        %--------
        begin = round(frame_length * (k - 1) / 2) + 1;
        endIndex =  begin + frame_length - 1;
        check = CheckVoice(voiceIndexArray,k); % Gọi hàm check vị trí là cùng tiếng nói
        %---------
        if check == 1 
            % Kiểm tra frame thứ k hiện tại có thuộc vùng tiếng nói hay
            % không
            range = begin : endIndex; %Chỉ số của các mẫu trong 1 frame
            frame = zeros(1, frame_length);
            for i = 1:frame_length % Duyệt chỉ mục frame
                if range(i) <= length(audioIn)
                    frame(i) =  hamm(i)*audioIn(range(i));
                end
            end
            P2 = abs(fft(frame, N)); %P2 gồm có 2 phần , một nửa ảo và một nửa thực 
            P_Voice = P2(1 : round(length(P2) / 2)); %Xét một nửa phía trái        
            freqV = linspace(0, Fs/2, length(P_Voice)); 
            %Tìm đỉnh bằng hàm findpeaks
            [y_value, y_peaks] =  findpeaks(P_Voice, freqV, 'NPeaks', 3, 'MinPeakDistance', 110, 'MinPeakHeight', 1);
            %---- Trường hợp 1: tìm được 3 peak ----%
            if length(y_peaks) == 3
                freq_max1 = abs(y_peaks(2) - y_peaks(1));     %Tìm hiệu giữa 2 cực đại cục bộ
                freq_max2 =abs(y_peaks(3) - y_peaks(2));
                if freq_max1 < 400 && freq_max1 > 70 && freq_max2 < 400 && freq_max2 > 70  %Kiểm tra tần số nằm trong khoảng (70-400)Hz
                    if freq_max1 >  1.5 * freq_max2
                        F0(index) = freq_max2;                 %Nếu freq_max1 lớn hơn 1.5 lần freq_max2 thì chỉ lấy freq_max2
                        F0_Voice = [F0_Voice freq_max2];
                    else
                        if freq_max2 >  1.5 * freq_max1     %Nếu freq_max2 lớn hơn 1.5 lần freq_max1 thì chỉ lấy freq_max1
                            F0(index) = freq_max1;
                            F0_Voice = [F0_Voice freq_max1];
                        else
                            F0(index) = (freq_max1+freq_max2)/2;% Lấy trung bình 2 tần số freq_max1 và freq_max2 suy ra F0
                            F0_Voice  = [F0_Voice (freq_max1+freq_max2)/2];
                        end
                    end
                end      
            else
                %---- Trường hợp 2: tìm được 2 peak thì lấy hiệu của 2 đỉnh là tần số cơ bản F0 ----%
                if length(y_peaks) == 2
                    freq_2peak = abs(y_peaks(1) - y_peaks(2));     %Tìm hiệu giữa 2 cực đại 
                    if freq_2peak<400 && freq_2peak>70                     %Kiểm tra tần số nằm trong khoảng (70-400)Hz
                        F0(index) = freq_2peak;                              %Neu co thi dua gia tri F0 do vao y_F0
                        F0_Voice = [F0_Voice freq_2peak];
                    end
                end
            end
            index = index + 1; 
        else
            if checkSilence == 1
                begin = round(frame_length * (k - 1) / 2) + 1;
                endIndex =  begin + frame_length - 1;               
                range = begin : endIndex; %Chỉ số của các mẫu trong 1 frame
                frame = zeros(1, frame_length);
                for i = 1:frame_length
                    if range(i) <= length(audioIn)
                        frame(i) =  hamm(i)*audioIn(range(i));
                    end
                end
                P2 = abs(fft(frame, N));%P2 gồm có 2 phần , một nửa ảo và một nửa thực 
                P_Unvoice = P2(1 : round(length(P2) / 2)); %Xét một nửa phía trái  
                freqU = linspace(1/Fs, Fs/2, length(P_Unvoice));
                checkSilence = 0;
            end        
           F0(index) = 0;
           index = index + 1;         
        end
    end
end
%--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------%

