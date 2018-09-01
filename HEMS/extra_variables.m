function [from_pv_demand,from_pv_battery,from_pv_grid,from_grid_battery,from_battery_demand,from_grid_demand,from_battery_grid] = extra_variables(NN,PV,PD,Pgplus,Pgminus,Pbplus,Pbminus,etaI,etaBc,etaBd)

Pres = etaI*PV - PD;
from_pv_demand = zeros(NN,1);
from_pv_battery = zeros(NN,1);
from_pv_grid = zeros(NN,1);
from_grid_battery = zeros(NN,1);
from_battery_demand = zeros(NN,1);
from_grid_demand = zeros(NN,1);
from_battery_grid = zeros(NN,1);
for i = 1:NN
    if PV(i) <= 0 %if PV is not generating
        from_pv_demand(i) = 0;
        from_pv_battery(i) = 0;
        from_pv_grid(i) = 0;
        
        if Pbplus(i) == 0
            if Pbminus(i) == 0
                from_battery_grid(i) = 0;
                from_battery_demand(i) = 0;
                from_grid_battery(i) = 0;
                from_grid_demand(i) = Pgplus(i) - from_grid_battery(i);
            elseif Pbminus(i) > 0
                if Pgplus(i) == PD(i)
                    from_grid_battery(i) = 0;
                    from_grid_demand(i) = Pgplus(i);
                    from_battery_demand(i) = 0;
                    from_battery_grid(i) = (1/etaBd)*etaI*Pbminus(i);
                elseif Pgplus(i) < PD(i)                  
                    if Pgminus(i) > 0
                        from_grid_battery(i) = 0;
                        from_battery_demand(i) = PD(i) - Pgplus(i);
                        from_battery_grid(i) = (1/etaBd)*etaI*Pbminus(i) - from_battery_demand(i);
                        from_grid_demand(i) = Pgplus(i);
                    else
                        from_battery_grid(i) = 0;
                        from_grid_battery(i) = 0;
                        from_battery_demand(i) = PD(i) - Pgplus(i); %OR from_battery_demand(i) = (1/etaBd)*etaI*Pbminus(i);
                        from_grid_demand(i) = Pgplus(i);
                    end
                else
                end
            end
            
        elseif Pbplus(i) > 0  %Pbminus(i) == 0
            from_grid_battery(i) = etaI*etaBc*Pbplus(i);
            from_grid_demand(i) = PD(i); 
            from_battery_grid(i) = 0;
            from_battery_demand(i) = 0;
        end
        
    else% If PV is generating
        if Pres(i) < 0  %But PV is less than demand
            if Pbplus(i) == 0 %if grid does not charge battery
                if Pbminus(i) == 0
                    from_pv_demand(i) = etaI*PV(i); %pv gives all it's got to demand
                    from_battery_grid(i) = 0;
                    from_battery_demand(i) = 0;
                    from_grid_battery(i) = 0;
                    from_grid_demand(i) = PD(i) - from_pv_demand(i);
                else
                    from_pv_demand(i) = etaI*PV(i); %pv gives all it's got to demand
                    from_pv_grid(i) = 0; % do not export to grid
                    from_pv_battery(i) = 0; %do charge the battery
                    from_grid_battery(i) = 0;
                    if Pgminus(i) == 0
                        from_battery_demand(i) = PD(i) - from_pv_demand(i); %battery brings the rest
                        from_battery_grid(i) = 0;
                    else
                        from_battery_demand(i) = PD(i) - from_pv_demand(i); %battery brings the rest
                        from_battery_grid(i) = Pgminus(i);
                    end
                end
            elseif Pbplus(i) > 0 %if grid is charging battery  Pbplus(i) > 0
                from_pv_demand(i) = etaI*PV(i);
                from_grid_demand(i) = PD(i) - from_pv_demand(i);
                from_grid_battery(i) = Pgplus(i) - from_grid_demand(i);
            end
        else % And PV is greater than demand
            from_pv_demand(i) = PD(i); %PV meets demand
            if Pbplus(i) == 0 %But if battery is full
                if Pbminus(i) == 0
                    from_battery_grid(i) = 0;
                    from_battery_demand(i) = 0;
                    from_pv_battery(i) = 0;
                    from_grid_battery(i) = 0;
                    from_pv_grid(i) = abs(Pres(i));
                else
                    from_grid_battery(i) = 0;
                    from_battery_demand(i) = 0;
                    from_pv_battery(i) = 0; %do not charge the battery
                    from_pv_grid(i) = abs(Pres(i)); %export to grid
                    from_battery_grid(i) = (1/etaBd)*etaI*Pbminus(i);
                end
            else %However, if battery is not full (Pbminus(i) == 0)
                from_pv_battery(i) = abs(Pres(i)); %charge the battery
                from_grid_battery(i) = etaI*etaBc*Pbplus(i) - from_pv_battery(i);
            end
            
        end
    end
end



end