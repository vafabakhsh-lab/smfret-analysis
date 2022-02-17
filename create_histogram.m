% Suppress no terminal ; error/warning
%#ok<*NOPTS>


close all;
fclose('all');
clear all;



path=input('where are the selected traces?  ', 's');
cd(path);
[FileName,PathName] = uigetfile('*.dat');
cd(PathName);

fret_bin=0.01; 
gap = -0.2:fret_bin:1.2;
temp_end=60;
N_smooth=4;

timeunit=input('Time unit: [default=0.03 sec] ');
    if isempty(timeunit)
        timeunit=0.03;
    end
    
% Older scripts did not keep track of the donor/acceptor background
% this variable will allow this script to be used to generate histograms
% from older data
number_columns = input(...
    'How many columns do the trace files have (4 or  6)?  ');
    if isempty(number_columns)
        number_columns=6;
    end

sg=max(size(gap));
N_His_f=zeros(1,sg);
N_His_fs=zeros(1,sg);
N_His_F=zeros(1,sg);
N_His_F_gamma=zeros(1,sg);
N_His_fret=zeros(1,sg);

h1=[];h2=[];h3=[];h4=[];h5=[];
N_total_XC=zeros(1,temp_end);
N_total_XC_old=zeros(1,temp_end);
N_total_XC_filter=zeros(1,temp_end);
N_total_XC_smooth=zeros(1,temp_end);
particle_count=0;
cr=[];u=1;

A=dir('*tr*.dat')
numberfiles=length(A) 

nm=0;

while (nm<numberfiles-0)
    nm=nm+1;
    A(nm).name;
    fname=A(nm).name(1:end-4);
    disp(fname)
    fid=fopen(A(nm).name,'r');
    if number_columns == 4
        dummy = fscanf(fid, '%g %g %g %g', [4 inf]);
    else
        dummy = fscanf(fid,'%g %g %g %g %g %g',[6 inf]);
    end
    fclose(fid);
    dec=size(dummy);
    points = dec(2);

    time=dummy(1,:)';
    donor=+dummy(2,:)';       
    acceptor=+dummy(3,:)';
    donorgamma=dummy(4,:)';
        

    
    I_f_don=[];I_b_don=[];I_f_acc=[];I_b_acc=[];C=[];f=[];b=[];
    I_f_fret=[];I_b_fret=[];
    smooth_donor=[];smooth_acceptor=[];
    fret=acceptor./(acceptor+donor);
    fret_smooth=smooth(fret,N_smooth);
    smooth_donor=smooth(donor,N_smooth);
    smooth_acceptor=smooth(acceptor,N_smooth);
        
    length=max(size(donor));
    n=[2];%2
    M=2;%2
    p=15;%15

    s=0;

    for N=n
        s=s+1;
        lower_bound=max(M,N+1);

        for i = lower_bound:length-N
            I_f_don(i,s)=(1/N)*sum(donor((i-N):i-1));
            I_b_don(i,s)=(1/N)*sum(donor(i+1:i+N));
    
            I_f_don_gamma(i,s)=(1/N)*sum(donorgamma((i-N):i-1));
            I_b_don_gamma(i,s)=(1/N)*sum(donorgamma(i+1:i+N));
    
            I_f_acc(i,s)=(1/N)*sum(acceptor((i-N):i-1));
            I_b_acc(i,s)=(1/N)*sum(acceptor(i+1:i+N));
    
            I_f_fret(i,s)=(1/N)*sum(fret((i-N):i-1));
            I_b_fret(i,s)=(1/N)*sum(fret(i+1:i+N));
        
        end


        for i = lower_bound:length-M-N+1
        f(i,s)=0;
        b(i,s)=0;
        fgamma(i,s)=0;
        bgamma(i,s)=0;
        ffret(i,s)=0;
        bfret(i,s)=0;
            for j=0:M-1
                j;
                f(i,s)=f(i,s)+((donor(i-j)-I_f_don((i-j),s))^2+(acceptor(i-j)-I_f_acc((i-j),s))^2);
                b(i,s)=b(i,s)+((donor(i+j)-I_b_don((i+j),s))^2+(acceptor(i+j)-I_b_acc((i+j),s))^2); 
        
                fgamma(i,s)=fgamma(i,s)+((donorgamma(i-j)-I_f_don_gamma((i-j),s))^2+(acceptor(i-j)-I_f_acc((i-j),s))^2);
                bgamma(i,s)=bgamma(i,s)+((donorgamma(i+j)-I_b_don_gamma((i+j),s))^2+(acceptor(i+j)-I_b_acc((i+j),s))^2);
        
                ffret(i,s)=ffret(i,s)+((fret(i-j)-I_f_fret((i-j),s))^2);  % filter the raw fret directly.
                bfret(i,s)=bfret(i,s)+((fret(i+j)-I_b_fret((i+j),s))^2);
            end
        f(i,s)=f(i,s)^(-p); b(i,s)=b(i,s)^(-p);
        fgamma(i,s)=fgamma(i,s)^(-p); bgamma(i,s)=bgamma(i,s)^(-p);
        ffret(i,s)=ffret(i,s)^(-p); bfret(i,s)=bfret(i,s)^(-p);

        end
    end

    C=(1./sum(f+b,2));     Cf=(1./sum(ffret+bfret,2));


    %%
    % $x^2+e^{\pi i}$
    f=bsxfun(@times,C,f);  ffret=bsxfun(@times,Cf,ffret);
    b=bsxfun(@times,C,b);  bfret=bsxfun(@times,Cf,bfret);

    Cgamma=(1./sum(fgamma+bgamma,2));     
    fgamma=bsxfun(@times,Cgamma,fgamma);  
    bgamma=bsxfun(@times,Cgamma,bgamma);  

    range=max(size(f));
    I_d_filter=0;I_a_filter=0;
    I_d_gamma_filter=0;
    fret_filter=0;

    for k=1:s
        I_d_filter= I_d_filter + f(1:range,k).* I_f_don(1:range,k) + b(1:range,k).*I_b_don(1:range,k);
        I_a_filter= I_a_filter + f(1:range,k).* I_f_acc(1:range,k) + b(1:range,k).*I_b_acc(1:range,k);
        I_d_gamma_filter= I_d_gamma_filter + fgamma(1:range,k).* I_f_don_gamma(1:range,k) + bgamma(1:range,k).*I_b_don_gamma(1:range,k);
    
        fret_filter= fret_filter + ffret(1:range,k).* I_f_fret(1:range,k) + bfret(1:range,k).*I_b_fret(1:range,k);
    end



    FRET_filter=I_a_filter./(I_a_filter+I_d_filter);
    FRET_filter_gamma=I_a_filter./(I_a_filter+I_d_gamma_filter); 

    lengthfret=max(size(fret));
    for t = 1:lengthfret-2
        if abs(fret(t+2)-fret(t+1))<0.02
            cr(u,1)=fret(t);
            cr(u,2)=fret(t+1);
            u=u+1;
        end
    end

    [n_f,~]=hist(fret,gap);          
    [n_fs,~]=hist(fret_smooth,gap);  
    [n_F,~]=hist(FRET_filter,gap);  
    [n_F_gamma,~]=hist(FRET_filter_gamma,gap);  
    [n_fret,x]=hist(fret_filter,gap);  
    % plot(x,n_F/sum(n_F),'r');hold on

        Ave_donor=mean(donor);        Ave_donor_filter=nanmean(I_d_filter);
        Ave_acceptor=mean(acceptor);  Ave_acceptor_filter=nanmean(I_a_filter);       
        std_donor=std(donor);         std_donor_filter=nanstd(I_d_filter);
        std_acceptor=std(acceptor);   std_acceptor_filter=nanstd(I_a_filter);
        delta_donor=donor-Ave_donor;  delta_donor_filter=I_d_filter-Ave_donor_filter;
        delta_acceptor=acceptor-Ave_acceptor; delta_acceptor_filter=I_a_filter-Ave_acceptor_filter;

        Ave_smooth_d = mean(smooth_donor); Ave_smooth_a = mean(smooth_acceptor);
        std_smooth_d= std(smooth_donor); std_smooth_a= std(smooth_acceptor);
        delta_smooth_d=smooth_donor - Ave_smooth_d; delta_smooth_a=smooth_acceptor - Ave_smooth_a;

        d_X_a=zeros(1,range+1);
        d_X_a_filter=zeros(1,range+1);
        d_X_a_smooth=zeros(1,range+1);

        delta_donor_filter(isnan(delta_donor_filter))=0; 
        delta_acceptor_filter(isnan(delta_acceptor_filter))=0;

        if range > temp_end
            for tau=0:range
            normalize_number=0;
            sum_tau=0;
            sum_tau_filter=0;
            sum_tau_smooth=0;
                for t=1:range-tau
                    sum_tau=sum_tau+(delta_donor(t)*delta_acceptor(t+tau));
                    sum_tau_filter=sum_tau_filter +(delta_donor_filter(t)*delta_acceptor_filter(t+tau));
                    sum_tau_smooth=sum_tau_smooth +(delta_smooth_d(t)*delta_smooth_a(t+tau));
                    normalize_number=normalize_number+1;
                end
                if normalize_number==0
                    d_X_a(tau+1)=NaN;
                else            
                    d_X_a_old(tau+1)=sum_tau/(normalize_number)/std_donor/std_acceptor;

                    d_X_a(tau+1)=sum_tau/(normalize_number)/(Ave_donor+Ave_acceptor);
                    d_X_a_smooth(tau+1)=sum_tau_smooth/(normalize_number)/(Ave_smooth_a+Ave_smooth_a);
                    d_X_a_filter(tau+1)=sum_tau_filter/(normalize_number)/(Ave_donor_filter+Ave_acceptor_filter);
                end
            end

            h1=subplot(5,2,[1,2,3,4]);

            cla(h1);cla(h2);cla(h3);cla(h4);    
            patchline(time,donor,'edgecolor','g','edgealpha',0.4);hold on;
            patchline(time(1:range),I_d_filter(1:range), 'linewidth',1.5,'edgecolor','g');hold on;
            patchline(time,acceptor,'edgecolor','r','edgealpha',0.4);hold on;
            patchline(time(1:range),I_a_filter(1:range), 'linewidth',1.5,'edgecolor','r');hold on;
            patchline(time,150+donor+acceptor,'edgecolor','k','edgealpha',0.4);hold on;
            patchline(time(1:range),150+I_a_filter(1:range)+I_d_filter(1:range), 'linewidth',1.5,'edgecolor','k');
            title(['File: ' fname '   Molecule:' int2str(nm)]);
            ylim([0 20+max(150+I_a_filter(1:range)+I_d_filter(1:range))]);
            xlim([min(time) max(time)]);
            grid on; zoom on;

            h2=subplot(5,2,[5,6]);
            patchline(time,fret,'edgecolor','c','edgealpha',0.5);hold on;
            patchline(time(1:range),FRET_filter(1:range),'linewidth',1.5,'edgecolor','b');hold on;
            ylim([-0.02 1.02]);
            xlim([min(time) max(time)]);
            grid on; zoom on;
            linkaxes([h1,h2], 'x');

            [~, indexAtMaxY] = max(n_f/sum(n_f));
            h3=subplot(5,2,[7,9]);
            plot(gap,n_f/sum(n_f),'*'); hold on;hold on;
            h4=subplot(5,2,[7,9]);
            plot(gap,N_His_f./sum(N_His_f),'r');
            xlim([-0.02 1.02]);
            legend(sprintf('Peak Center = %0.3f',gap(indexAtMaxY(1))));

            h5=subplot(5,2,[8,10]);
            cla(h5);
            plot(d_X_a_old(1:temp_end),'*b');hold on;
            plot(N_total_XC_old/particle_count,'r');
            ylim([-0.6 0.2]);
            legend(sprintf('XC magnitude = %0.3f',d_X_a_old(2)));


            option=input('include (enter), skip (s) ', 's');
            if option=='s'
                N_His_f=N_His_f;
                N_His_fs=N_His_fs;
                N_His_F=N_His_F;
                N_His_F_gamma=N_His_F_gamma;
                N_His_fret=N_His_fret;
            else
                N_His_f=N_His_f+n_f./sum(n_f);
                N_His_fs=N_His_fs+n_fs./sum(n_fs);
                N_His_F=N_His_F+n_F./sum(n_F); 
                if donorgamma(1)>0 
                    N_His_F_gamma=N_His_F_gamma+n_F_gamma./sum(n_F_gamma); 
                end
                N_His_fret=N_His_fret+n_fret./sum(n_fret);
        %%%%%%%%%%%%%%% CORRELATION  %%%%%%%%%%%%%%%%

                N_total_XC=N_total_XC + d_X_a(1:temp_end);
                N_total_XC_filter=N_total_XC_filter + d_X_a_filter(1:temp_end);
                N_total_XC_smooth=N_total_XC_smooth + d_X_a_smooth(1:temp_end);
                N_total_XC_old=N_total_XC_old + d_X_a_old(1:temp_end);
                
                
                particle_count=particle_count+1;


            end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        end
end

% Normalize data by particle counts
N_total_XC = N_total_XC/particle_count;
N_total_XC_old = N_total_XC_old/particle_count;
N_total_XC_filter = N_total_XC_filter/particle_count;

N_His_f = N_His_f/particle_count;
N_His_fs = N_His_fs/particle_count;
N_His_fret = N_His_fret/particle_count;
N_His_F = N_His_F/particle_count;
N_His_F_gamma = N_His_F_gamma/particle_count;

%%%%%        TIME                     No filter                  old normalization              filtered data  
output=[timeunit*(0:temp_end-1)' (N_total_XC)' (N_total_XC_old)' (N_total_XC_filter)'];
time_fname=['total_XC.dat'];
save(time_fname,'output','-ascii') ;

%xcHeaders = {"TIME", "NoFilter", "OldNorm", "FilteredData"};
%csvXC = vertcat(xcHeaders, num2cell(output));
csvwrite('total_XC.csv', output);

%%%%% TIME  raw data  smoothed fret  filtered fret   filtered_donor/filtered _acceptor 
output=[x' , N_His_f' N_His_fs' N_His_fret' N_His_F' N_His_F_gamma'];
time_fname=['total_FRET.dat'];
save(time_fname,'output','-ascii') ;

%fretHeaders = {"FRET", "Raw", "Smoothed", "Filtered", "FilteredD.A", "gamma"} ;
%csvFRET = vertcat(fretHeaders, num2cell(output));
csvwrite('total_FRET.csv', output);








%%%% ChangeLog %%%%
% 4/2/2019 Michael Schamber
% Changed plotting windows to give extra panel to XC and changed limits so
% the XC plot is more visible. This required making the FRET trace plot a
% little thinner so that it fits. The FRET histogram plot remains the same
% size.

