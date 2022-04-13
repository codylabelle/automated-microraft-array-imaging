global user sys
T=1:1:24;
for i = 1:1:24
    VM(i) = [sys{i}.VirtualAddressSpace.Available];
    PM(i) = [sys{i}.PhysicalMemory.Available/(10^9)];
    UM(i) = [user{i}.MemUsedMATLAB/(10^9)];
    VMT(i) = [sys{i}.VirtualAddressSpace.Total/(10^9)];
end
%Plot virtual memory available
plot(T,VM)
ylim([0 1000000]);

%Plot virtual memory available
plot(T,VMT)
ylim([0 1000000000000000]);

%Plot Physical memory available
plot(T,PM)
ylim([0 50]);

%Plot Memory used by matlab
plot(T,UM)
ylim([0 50]);

