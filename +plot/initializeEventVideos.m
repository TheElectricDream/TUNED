function [hFigs, hAxs, hImgs, videoWriters] = initializeEventVideos(cohOut, atsOut, imgSz)

    if atsOut

        % Set up the ATS output figure
        hFigATS = figure('Visible','off','Position',[100 100 imgSz]);  
        hAxATS  = axes('Parent',hFigATS);                                   
        colormap(hAxATS,'gray');                                          
        colorbar(hAxATS);                                              
        set(hAxATS,'FontSize',16,'Color','white');
        set(hAxATS, 'XTick', [], 'YTick', []);
        initial_data = nan(imgSz(2), imgSz(1)); 
        hImgATS = imagesc(hAxATS, initial_data); 
        set(hImgATS, 'AlphaData', ~isnan(initial_data));
        
        % Set up video writer
        videoFileName = '/home/alexandercrain/Videos/Research/adaptive_time_surface_accumulator.avi';
        videoWriters{1} = VideoWriter(videoFileName);
        videoWriters{1}.FrameRate = 60;  % Set the frame rate
        open(videoWriters{1});

    else
        
        % Set handles to be empty
        hFigATS = [];
        hAxATS = [];
        hImgATS = [];

    end

    if cohOut

        % Set up the ATS output figure
        hFigCOH = figure('Visible','off','Position',[100 100 imgSz]);  
        hAxCOH  = axes('Parent',hFigCOH);                                   
        colormap(hAxCOH,'gray');                                          
        colorbar(hAxCOH);                                              
        set(hAxCOH,'FontSize',16,'Color','white');
        initial_data = nan(imgSz(2), imgSz(1)); 
        hImgCOH = imagesc(hAxCOH, initial_data); 
        set(hImgCOH, 'AlphaData', ~isnan(initial_data));

        % Set up video writer
        videoFileName = '/home/alexandercrain/Videos/Research/coherence_map_output.avi';
        videoWriters{2} = VideoWriter(videoFileName);
        videoWriters{2}.FrameRate = 60;  % Set the frame rate
        open(videoWriters{2});

    else
        
        % Set handles to be empty
        hFigCOH = [];
        hAxCOH = [];
        hImgCOH = [];

    end


    % Store handles in output variables
    hFigs{1} = hFigATS; 
    hAxs{1} = hAxATS; 
    hImgs{1} = hImgATS;

    % Store handles for coherence output
    hFigs{2} = hFigCOH; 
    hAxs{2} = hAxCOH; 
    hImgs{2} = hImgCOH;

end
