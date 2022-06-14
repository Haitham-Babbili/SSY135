%% Time and Frequency-Varying Rician Fading Channels
clc , clear, close all
Ts = 0.1e-3;                    %symbol time Ts in s
Ns = 1e6;                       %number of samples per channel
N = 300;
M = 64;
fDTs =[0.1,0.005];
L = [1 2 3];
kc =[0 1 10];

for k=1:3 % For all kc
    for j=1:2 % For all fDTs
        %L = 1
        c_1(1,:) = channelByFilter(Ts,Ns,N,fDTs(j)/Ts) +kc(k)*ones(1,1000000); % channel by filter + kc
        c_1(2,:) = zeros(1,1000000); % Can not make it work without this vector in c_1
        c_1_padded=[c_1 zeros(2,M-L(1))]; % Pad with M-L zeros
        C_1 = fft(c_1_padded,M);
        C_1 = C_1(1:M,1:N);
        
        %L = 2
        c_2 = zeros(L(2),Ns);
        for l=1:L(2)
            temp = channelByFilter(Ts,Ns,N,fDTs(j)/Ts) +kc(k)*ones(1,1000000); % channel by filter + kc
            c_2(l,:) = temp;
        end
        c_2_padded = [c_2 zeros(2,M-L(2))]; % Pad with M-L zeros
        C_2 = fft(c_2_padded,M);
        C_2 = C_2(1:M,1:N);
        
        % %L = 3
        c_3 = zeros(L(3),Ns);
        for l=1:L(3)
            temp = channelByFilter(Ts,Ns,N,fDTs(j)/Ts) +kc(k)*ones(1,1000000); % channel by filter + kc
            c_3(l,:) = temp;
        end
        c_3_padded=[c_3 zeros(3,M-L(3))]; % Pad with M-L zeros
        C_3 = fft(c_3_padded,M);
        C_3 = C_3(1:M,1:N);
        
        figure('Renderer', 'painters', 'Position', [10 10 1600 1000])
        colormap('jet')
      %%% sgtitle(['k = ' + string(kc(k))+ ', fDTs = ' + string(fDTs(j))])
        %plot for L=1
        subplot(3,2,1)
        mesh(0:N-1,0:M-1, abs(C_1));  xlabel('time[T_s]'), ylabel('frequency[1/NT_s]|'), zlabel('|C[k,n]|')
        subplot(3,2,2)
        surf(1:N, 0:M-1,abs(C_1), 'MeshStyle', 'row'); xlabel('time[T_s]'), ylabel('frequency[1/NT_s]|'), view(0,90)
        %plot for L=2
        subplot(3,2,3)
        mesh(0:N-1,0:M-1, abs(C_2));  xlabel('time[T_s]'), ylabel('frequency[1/NT_s]|'), zlabel('|C[k,n]|')
        subplot(3,2,4)
        surf(1:N, 0:M-1,abs(C_2), 'MeshStyle', 'row'); xlabel('time[T_s]'), ylabel('frequency[1/NT_s]|'), view(0,90)
        %plot for L=3
        subplot(3,2,5)
        mesh(0:N-1,0:M-1, abs(C_3));  xlabel('time[T_s]'), ylabel('frequency[1/NT_s]|'), zlabel('|C[k,n]|')
        subplot(3,2,6)
        surf(1:N, 0:M-1,abs(C_3), 'MeshStyle', 'row'); xlabel('time[T_s]'), ylabel('frequency[1/NT_s]|'), view(0,90)
    end
end
