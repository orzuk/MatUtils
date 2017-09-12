% Solve Einstein's riddle

colors = {'red', 'green', 'blue', 'yellow', 'white'};
nations = {'england', 'denmark', 'norway', 'sweden', 'germany'};
cigars = {'pall mall', 'dunhills', 'blend', 'prince', 'blue masters'};
pets = {'birds', 'dogs', 'horses', 'cats', 'fish'};
drinks = {'tea', 'milk', 'water', 'coffee', 'bier'};

P = perms(1:5);
n = length(P);

good_color_inds = zeros(n,1);
for i_color = 1:n
    permed_colors = colors(P(i_color,:));
    % 1. check if green house is just left (smaller by 1) than white house
    green_ind = strmatch('green', permed_colors);
    white_ind = strmatch('white', permed_colors);
    if(white_ind == green_ind+1)
        good_color_inds(i_color)=1;
    end
end
good_color_inds = find(good_color_inds)

good_nation_inds = zeros(n,1);
for i_nation = 1:n
    permed_nations = nations(P(i_nation,:));
    % 2. check if norway is first house
    if(strmatch('norway', permed_nations) == 1)
        good_nation_inds(i_nation)=1;
    end
end
good_nation_inds = find(good_nation_inds)

for i_color = good_color_inds'
    permed_colors = colors(P(i_color,:));
    for i_nations = good_nation_inds'
        permed_nations = nations(P(i_nations,:));
        if(abs( strmatch('norway', permed_nations) - ...
                strmatch('blue', permed_colors) )~=1 ) % 3. check norway next to blue
            continue;
        end
        if(strmatch('england', permed_nations) ~= ...
                strmatch('red', permed_colors) ) % 4. check english in red
            continue;
        end
        
        
        %   sprintf(' Run:')
        %   good_nations = permed_nations
        %   good_colors = permed_colors
        
        
        for i_cigars = 1:n
            permed_cigars = cigars(P(i_cigars,:));
            if(strmatch('yellow', permed_colors) ~= ...
                    strmatch('dunhills', permed_cigars) ) % 5. dunhills in yellow
                continue;
            end
            
            if(strmatch('germany', permed_nations) ~= ...
                    strmatch('prince', permed_cigars) ) % 6. german smokes prince
                continue;
            end
            
            for i_pets = 1:n
                permed_pets = pets(P(i_pets,:));
                if(strmatch('sweden', permed_nations) ~= ...
                        strmatch('dogs', permed_pets) ) % 7. swede keeps dogs
                    continue;
                end
                if(strmatch('birds', permed_pets) ~= ...
                        strmatch('pall mall', permed_cigars) ) % 8. pall mall smoker keeps birds
                    continue;
                end
                
                if(abs( strmatch('blend', permed_cigars) - ...
                        strmatch('cats', permed_pets) )~=1 ) % 9. blend smoker next to cats
                    continue;
                end
                if(abs( strmatch('dunhill', permed_cigars) - ...
                        strmatch('horses', permed_pets) )~=1 ) % 10. dunhill smoker next to horses
                    continue;
                end
                for i_drinks = 1:n
                    permed_drinks = drinks(P(i_drinks,:));
                    % Check constraints
                    if(strmatch('denmark', permed_nations) ~= ...
                            strmatch('tea', permed_drinks) ) % 11. dane drinks tea
                        continue;
                    end
                    if(strmatch('green', permed_colors) ~= ...
                            strmatch('coffee', permed_drinks) ) % 12. green drinks coffee
                        continue;
                    end
                    if(strmatch('milk', permed_drinks) ~= 3) % 13. center house drinks milk
                        continue;
                    end
                    if(strmatch('blue master', permed_cigars) ~= ...
                            strmatch('bier', permed_drinks) ) % 14. blue-master drinks bier
                        continue;
                    end
                    if(abs( strmatch('blend', permed_cigars) - ...
                            strmatch('water', permed_drinks) )~=1 ) % 10. dunhill smoker next to horses
                        continue;
                    end
                    
                    sprintf('Found Good Assignment!!!')
                    colors_are = permed_colors
                    nations_are = permed_nations
                    pets_are = permed_pets
                    cigars_are = permed_cigars
                    drinks_are = permed_drinks
                    
                    
                end
            end
        end
    end
end

