@testset "TestInputOutput: Test internal consistency of input/output condition calculations" begin

    @testset "Csys" begin
        
        @testset "Temperature Effect" begin
            c1 = whole_system(pHtot=8.1, TA=2300, T_in=20, T_out=30)
            c2 = whole_system(pHtot=c1.pHtot, TA=2300, T_in=30, T_out=20)
            @test c1.pHtot_in ≈ c2.pHtot atol=1e-6
        end
        
        @testset "Salinity Effect" begin
            c1 = whole_system(pHtot=8.1, TA=2300, S_in=28.2, S_out=38.1)
            c2 = whole_system(pHtot=c1.pHtot, TA=2300, S_in=38.1, S_out=28.2)
            @test c1.pHtot_in ≈ c2.pHtot atol=1e-6
        end
        
        @testset "Pressure Effect" begin
            c1 = whole_system(pHtot=8.1, TA=2300, P_in=0, P_out=400)
            c2 = whole_system(pHtot=c1.pHtot, TA=2300, P_in=400, P_out=0)
            @test c1.pHtot_in ≈ c2.pHtot atol=1e-6
        end

    end

    @testset "CBsys" begin
        
        @testset "Temperature Effect" begin
            c1 = whole_system(pHtot=8.1, TA=2300, T_in=20, T_out=30, δBT=39.4)
            c2 = whole_system(pHtot=c1.pHtot, TA=2300, T_in=30, T_out=20, δBT=39.4)
            @test c1.pHtot_in ≈ c2.pHtot atol=1e-6
        end
        
        @testset "Salinity Effect" begin
            c1 = whole_system(pHtot=8.1, TA=2300, S_in=28.2, S_out=38.1, δBT=39.4)
            c2 = whole_system(pHtot=c1.pHtot, TA=2300, S_in=38.1, S_out=28.2, δBT=39.4)
            @test c1.pHtot_in ≈ c2.pHtot atol=1e-6
        end
        
        @testset "Pressure Effect" begin
            c1 = whole_system(pHtot=8.1, TA=2300, P_in=0, P_out=400, δBT=39.4)
            c2 = whole_system(pHtot=c1.pHtot, TA=2300, P_in=400, P_out=0,δBT=39.4)
            @test c1.pHtot_in ≈ c2.pHtot atol=1e-6
        end

    end

end