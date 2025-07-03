Idx = [];
for i = 1:5298
    Addmissiblet = TrainConfiguration(i).Addmissiblet;
    if Addmissiblet(1) < 36000 && Addmissiblet(2) > 36900
        Idx = [Idx i];
        if length(Idx) >= 100
            break
        end
    end
end
TrainConfiguration_example = TrainConfiguration(Idx);
TrainConfiguration_example = rmfield(TrainConfiguration_example, 'xKvturn');
TrainConfiguration_example = rmfield(TrainConfiguration_example, 'ILO');
save('TrainConfiguration_example.mat','TrainConfiguration_example');