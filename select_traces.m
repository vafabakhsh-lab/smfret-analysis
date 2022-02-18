function trace_selection
%by Hajin to build FRET histogram of selected range of traces with background correction 1/7/2010
%updated to be able to 1) delete the most recent selection, 2) export
%raw data, and 3) go back one molecule 1/20/2010

% now modified by reza to get specific file and changed the controls
% this is the las version. It averages all the FRET values in the selected
% segment and gives one number. each segment contributes one data point.

% final modification: now saves  the average and full segment data and also sm average in
% different files. 
% each file has 2 columns: total intensity and fret
% changed the averaging of histograms


% Modified by Reza March 31 to correct for Gamma

% Excess code removed and variables renamed by Michael Schamber on Feb 2022
% additional comments made to try and explain code more clearly

% Suppress no terminal ; error/warning
%#ok<*NOPRT>

close all;
fclose('all');
clear all;
% request file directory of data to be analyzed
pth=input('Directory: ', 's'); % 's' - autoconverts input to string
cd(pth); % directory so saved traces are in same location


sm=3; % specifies number of data points to average
gamma = 1.0; % specifies the default Gamma factor
bin_size = 0.01; % specifies the size of fret bins
fret_bins = -0.2 : bin_size : 1.2; % specifies range for fret bins
sg=size(fret_bins);
sg=sg(2);
N_His=zeros(1,sg);
N_His_raw=zeros(1,sg);

fname=input('index # of filename [default=1]  '); 
if isempty(fname)
fname=1;
end
fname=num2str(fname);
%filename=['film' fname '.traces'] 
fid=fopen(['film' fname '.traces'],'r');
 
A=dir;
[nf,dumn]=size(A);
numberoffiles=length(A)-2;
disp(['Number of files: ',num2str(numberoffiles)]);

timeunit=input('Time unit: [default=0.03 sec] ');
    if isempty(timeunit)
        timeunit=0.03;
    end
    
    
leakage=input('Donor leakage correction: [default=0.085] ');
if isempty(leakage)
  leakage=0.085;
end
    
    
bg_definition_input=input('Manually define background? y or n [default=y] ');
if isempty(bg_definition_input)
    bg_definition_input='y';
end

bg_definition=1;

if bg_definition_input=='n'
    bg_definition=0;
end


number_of_frames=0;
total_N_of_frames=0;
was_deleted=0;


% repeat over all the files in the directory



    
len=fread(fid,1,'int32'); 
disp(['The len of the time traces is ',num2str(len)]);
Ntraces=fread(fid,1,'int16');
disp(['The number of traces is ',num2str(Ntraces/2)]);

time=(0:(len-1))*timeunit;

raw=fread(fid,Ntraces*len,'int16');
fclose(fid);
% convert into donor and acceptor traces
index=(1:Ntraces*len);
Data=zeros(Ntraces,len);
donor=zeros(Ntraces/2,len);
acceptor=zeros(Ntraces/2,len);
Data(index)=raw(index);

for j=1:(Ntraces/2)
    donor(j,:)=Data(j*2-1,:);
    acceptor(j,:)=Data(j*2,:);
end

% define FRET
fret=zeros(Ntraces/2,len);
for j=1:Ntraces/2
    for k=1:len
        if donor(j,k)+acceptor(j,k)==0
            fret(j,k)=0.5;
        else
            fret(j,k)=gamma.*(acceptor(j,k)-leakage*donor(j,k))/(donor(j,k)+gamma.*acceptor(j,k));
        end
        if fret(j,k)>1.5
            fret(j,k)=1.5;
        end
        if fret(j,k)<-0.5
            fret(j,k)=-0.5;
        end
    end
end

hdl=gcf;        
set(hdl,'Position',[4,300,1276,648]);
    
%%%%% JUMP TO SPECIFIC MOLECULE
%%%%% j=250;

j=0;

while(j<Ntraces/2)
    j=j+1;

    clf % clear figure to prevent lag caused by plotting

    % plot donor, acceptor, and FRET
    h1=subplot(5,1,[1,2,3]);
    plot(time,donor(j,:),'g',time,(acceptor(j,:))-abs(leakage*donor(j,:)),'r', time,(acceptor(j,:)+donor(j,:)+200), 'k' );
    ylim([-25 max((acceptor(j,:)+donor(j,:)+400))]);
    xlim([-20*timeunit max(time)+100*timeunit]);

    title(['File:' fname '   Molecule:' int2str(j-1)]);
    grid on;
    zoom on;
    h2=subplot(2,1,2);
    plot(time,fret(j,:),'b');
     ylim([-0.02 1.02]);
     xlim([-20*timeunit max(time)+100*timeunit]);

    grid on;
    zoom on;
    linkaxes([h1,h2], 'x');
    %frames_to_avg=3;

    option=input('skip (enter), Work on this trace(s), go back (b), save trace (t) ','s');
    if option=='t'
        fname1=[fname ' tr' num2str(j-1) '.dat'];
        output=[time' donor(j,:)' acceptor(j,:)'];
        save(fname1,'output','-ascii') ;
        j=j-1;
    else

%         if option=='s'
%             disp('Molecule skipped');
%         else
        if option=='b'
            if j==1
                disp('This is the first molecule');
                j=j-1;
            else
                disp('Going to previous molecule');
                j=j-2;
            end
        else
            if option=='c'
                if (j==1 && i==3)
                    disp('This is the first molecule');
                else
                    if number_of_frames==0
                        disp('There is no selection to delete');
                    else
                        if was_deleted==0
                            total_N_of_frames=total_N_of_frames-number_of_frames;                        
                            disp('Previous selection was cancelled');
                            was_deleted=1;
                        else
                            disp('Previous selection was already cancelled');
                        end
                    end
                end
                j=j-1;
            else
                if option=='s'
                disp('Click twice to select the range (press enter to skip)');
                [x,~]=ginput(2);

                if isempty(x) || size(x)==1
                    x=ones(2);
                else
                    x(1)=round(x(1)/timeunit);
                    x(2)=round(x(2)/timeunit);
                    if x(1)<1
                        x(1)=1;
                    end
                    if x(2)<1
                        x(2)=1;
                    end
                    if x(1)>len
                        x(1)=len;
                    end
                    if x(2)>len
                        x(2)=len;
                    end
                end

                if x(1)>x(2)
                    temp=x(1);
                    x(1)=x(2);
                    x(2)=temp;
                end

                disp(['Selected range is frame ',num2str(x(1)),' to frame ',num2str(x(2))]);
                number_of_frames=x(2)-x(1)+1;
                was_deleted=0;


                x(2)=x(2)-mod((x(2)-x(1)),sm);



                % Manual background selection

                if bg_definition
                    disp('Click 4X  to select the range to calculate background (enter to skip)');
                    [xbc3,~]=ginput(2);
                    [xbc5,~]=ginput(2);

                    if xbc3(2)> xbc3(1)
                        xbc3(1)=round(xbc3(1)/timeunit);
                        xbc3(2)=round(xbc3(2)/timeunit);
                        bg_donor=mean(donor(j,xbc3(1):xbc3(2)));
                    else
                        bg_donor=input('donor background? ');
                    end

                    if xbc5(2)> xbc5(1)
                        xbc5(1)=round(xbc5(1)/timeunit);
                        xbc5(2)=round(xbc5(2)/timeunit);
                        bg_acceptor=mean(acceptor(j,xbc5(1):xbc5(2)) - leakage.*donor(j,xbc5(1):xbc5(2)));
                    else
                        bg_acceptor=input('acceptor background? ');
                    end




                end

                disp(['Backgrounds of donor and acceptor are ',num2str(bg_donor),' and ',num2str(bg_acceptor)]);

               %%%% GAMMA FACTOR CORRECTION  %%%%%%%%%%%%%
               %%%%%% Either manually enter GAMMA or click 4 times (donor low, donor high then acceptor low, acceptor high)
               %%%%%%%%%%% let the program measure Gammma. If you don't know Gamma put 0.%%%

               gamma=0;
               deltaA=0;
               deltaD=0; 

                 fid = fopen(['gamma' fname '.dat'],'a+');
                 fprintf(fid,'%4.3f\n', gamma);
                 fclose(fid);
                %%%% END of GAMMA FACTOR CORRECTION%%%%
                %%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                if number_of_frames>1


                    temp_donorsm=smooth(donor(j,x(1):x(2)),3)-bg_donor;
                    temp_acceptorsm=smooth(acceptor(j,x(1):x(2)),3)-bg_acceptor;
                    allfretsm=(temp_acceptorsm-abs(leakage.*temp_donorsm))./(temp_acceptorsm+temp_donorsm);
                    [nn,xout] = hist(allfretsm,fret_bins);
                    nn_norm=nn/sum(nn);
                    N_His=(N_His+nn_norm)/sum(N_His+nn_norm);
                    N_His_raw=N_His_raw+nn;
                    
                    %%%%%%%%%%%%%%
                    output=[xout' , N_His', N_His_raw'];
                    time_fname=['selectedFRET_SMave' fname '.dat'];
                    save(time_fname,'output','-ascii') ;



                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    mean_donor=mean(donor(j,x(1):x(2)))-bg_donor;
                    mean_acceptor=mean(acceptor(j,x(1):x(2)))-bg_acceptor;
                    allfret=(mean_acceptor-abs(leakage*mean_donor))/(mean_acceptor+mean_donor)
                    allfret_gamma=(mean_acceptor-abs(leakage*mean_donor))/(mean_acceptor+gamma*mean_donor);
                    selfret=(acceptor(j,x(1):x(2))-bg_acceptor-abs(leakage*(donor(j,x(1):x(2)))))./((donor(j,x(1):x(2)))+acceptor(j,x(1):x(2))-bg_donor-bg_acceptor);
                    seltotal=donor(j,x(1):x(2))+acceptor(j,x(1):x(2));
                    totalI=(mean_acceptor+mean_donor);

                    fid = fopen(['selectedFRETavseg' fname '.dat'],'a+');
                    fprintf(fid, '%4.0f\t %5.1f\t  %4.3f\t %4.3f\n', j,totalI,allfret,gamma);
                    fclose(fid);
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    selres=[seltotal ; selfret];
                    fid = fopen(['selectedFRET' fname '.dat'],'a+');
                    fprintf(fid,'%5.1f  %4.3f\n', selres);
                    fclose(fid);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%This is the selected traces that we will use for future analysis. 
%%% the columns are: TIME, DONOR (background corrected), Acceptor 
%%% (background and leakage corrected) and gamma*DONOR. In case of
%%% Gamma, FRET=A/(A+gamma*D).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                    ttsize=max(size(time(x(1):x(2))));
                    fname5=[fname ' tr' num2str(j) '.dat'];
                    output5=[time(x(1):x(2))'  (donor(j,x(1):x(2))-bg_donor)' (acceptor(j,x(1):x(2))-abs(leakage*(donor(j,x(1):x(2))))-bg_acceptor)' gamma*(donor(j,x(1):x(2))-bg_donor)' bg_donor*ones(ttsize,1) bg_acceptor*ones(ttsize,1)];
                    save(fname5,'output5','-ascii') ;


                    total_N_of_frames=total_N_of_frames+number_of_frames;
                    allfret=[];

                    more_range=input('Press b to select more range [default=no] ','s');
                    if more_range=='b'
                        j=j-1;
                    end
                end

                end
            end
        end
    end
end
    


