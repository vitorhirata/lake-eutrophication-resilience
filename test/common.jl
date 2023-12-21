@testset "#_entropy" begin
    @testset "Low P0 and number of decisions one" begin
        P0 = 0.5
        I = 0.1
        decision_step = 5.0
        number_decision = 1
        deterministic = true
        maximum_options = 10

        entropy = PathwayDiversity._entropy(P0, I, decision_step, number_decision, deterministic, maximum_options)
        number_options = PathwayDiversity.number_possible_influx(P0, maximum_options)

        @test number_options == 10
        @test entropy ≈ 2.302585 rtol=1e-5
    end

    @testset "High P0 and number of decisions one" begin
        P0 = 2.0
        I = 0.1
        decision_step = 5.0
        number_decision = 1
        deterministic = true
        maximum_options = 10

        entropy = PathwayDiversity._entropy(P0, I, decision_step, number_decision, deterministic, maximum_options)
        number_options = PathwayDiversity.number_possible_influx(P0, maximum_options)

        @test number_options == 4
        @test entropy ≈ 1.386294 rtol=1e-5
    end

    @testset "Low P0 and number of decisions two" begin
        P0 = 0.5
        I = 0.1
        decision_step = 5.0
        number_decision = 2
        deterministic = true
        maximum_options = 10

        entropy = PathwayDiversity._entropy(P0, I, decision_step, number_decision, deterministic, maximum_options)

        @test entropy ≈ 4.416255 rtol=1e-5
    end

    @testset "High P0 and number of decisions two" begin
        P0 = 2.0
        I = 0.1
        decision_step = 5.0
        number_decision = 2
        deterministic = true
        maximum_options = 10

        entropy = PathwayDiversity._entropy(P0, I, decision_step, number_decision, deterministic, maximum_options)

        @test entropy ≈ 2.756454 rtol=1e-5
    end

    @testset "High P0, number of decisions five and initial influx high" begin
        P0 = 0.5
        I = 0.15
        decision_step = 5.0
        number_decision = 6
        deterministic = true
        maximum_options = 10

        entropy = PathwayDiversity._entropy(P0, I, decision_step, number_decision, deterministic, maximum_options)

        @test entropy ≈ 11.856449 rtol=1e-5
    end
end
