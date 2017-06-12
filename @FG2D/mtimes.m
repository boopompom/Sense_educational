function res = mtimes(fg,data) 
% Fessler 2D NUFFT operator working on 12D reconframe data
% Data and kspace come in as cells, which represent different data chunks
% Try indexing in parfor loop only on first entry, probably saves speed
%
% Tom Bruijnen - University Medical Center Utrecht - 201704 




    % Check what dimensions require new trajectory coordinates
    Id=fg.Id;
    Kd=fg.Kd;

    if fg.adjoint==1    % non-Cartesian k-space to Cartesian image domain || type 1
        
        % Reshape data that goes together into the nufft operator
        data=reshape(data,[Kd(1)*Kd(2) Kd(3)]);

        % Preallocate response cell
        res=[];

        % Track progress
        if fg.verbose;parfor_progress(prod(Kd(3:end)));end
        
        % Loop over all dimensions and update k if required
        % For now I assumed that different Z always has the same trajectory

            % Convert data to doubles, required for the function
            data_tmp=double(data);

            % Parallize over the receivers (always has same traj)
            for coil=1:Kd(3)
                % Save in temporarily matrix, saves indexing time
                res_tmp(:,coil)=matrix_to_vec(nufft_adj(data_tmp(:,coil),...
                    fg.st)/sqrt(prod(fg.Id(1:2))));      
                
                % Track progrss
                if fg.verbose;parfor_progress;end
            end
            
            % Store output from all receivers
            res=reshape(res_tmp,Id);
            
        % Reset progress file
        if fg.verbose;parfor_progress(0);end

    else         % Cartesian image domain to non-Cartesian k-space || type 2

        % Preallocate response cell
        res=[];
        
        % Loop over all dimensions and update k if required
        % For now I assumed that different Z always has the same trajectory

                % Convert data to doubles, required for the function
                data_tmp=double(data);
                
                % Parallize over the receivers (always has same traj)
                for coil=1:Kd(3)
                    % Save in temporarily matrix, saves indexing time
                    res_tmp(:,coil)=matrix_to_vec(nufft(data_tmp(:,:,coil),...
                        fg.st)/sqrt(prod(fg.Id(1:2))));
                end

                % Store output from all receivers
                res=reshape(res_tmp,Kd);

end

% END  
end 
