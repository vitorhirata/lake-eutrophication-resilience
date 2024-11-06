@testset "#_entropy standard calculus" begin
    @testset "Low P0 and number of decisions one" begin
        P0 = 0.5
        I = 0.1
        decision_step = 5.0
        number_decision = 1
        deterministic = true
        max_options = 10

        entropy = PathwayDiversity._entropy(P0, I, decision_step, number_decision;
                                            deterministic=deterministic, max_options=max_options)
        number_options = PathwayDiversity.number_possible_influx(P0, max_options)

        # Giving the assumption all decision have the same probability, we can compute the entropy using:
        # entropy = - number_options * (1/number_options)* log(1/number_options) = log(number_options)
        @test number_options == 10
        @test entropy ≈ log(number_options) rtol=1e-5
    end

    @testset "High P0 and number of decisions one" begin
        P0 = 2.0
        I = 0.1
        decision_step = 5.0
        number_decision = 1
        deterministic = true
        max_options = 10

        entropy = PathwayDiversity._entropy(P0, I, decision_step, number_decision;
                                            deterministic=deterministic, max_options=max_options)
        number_options = PathwayDiversity.number_possible_influx(P0, max_options)

        @test number_options == 4
        @test entropy ≈ log(number_options) rtol=1e-5
    end

    @testset "Low P0 and number of decisions two" begin
        P0 = 0.5
        I = 0.1
        decision_step = 5.0
        number_decision = 2
        deterministic = true
        max_options = 10

        entropy = PathwayDiversity._entropy(P0, I, decision_step, number_decision;
                                            deterministic=deterministic, max_options=max_options)

        # To compute the entropy in this case the number of influx from the last state will simplify with the
        # probability of choosing the influx. So, from each final state we will have:
        # last_n_opt*(1/first_n_opt*last_n_opt)*log(1/first_n_opt*last_n_opt)
        # Ex. -9*(1/9*10)*log(1/90) = (1/10)*log(90)
        @test entropy ≈ (5*log(100)+log(90)+log(80)+log(70)+log(60)+log(50))/10 rtol=1e-5
    end

    @testset "High P0 and number of decisions two" begin
        P0 = 2.0
        I = 0.1
        decision_step = 5.0
        number_decision = 2
        deterministic = true
        max_options = 10

        entropy = PathwayDiversity._entropy(P0, I, decision_step, number_decision;
                                            deterministic=deterministic, max_options=max_options)

        @test entropy ≈ (log(5*4)+log(4*4)+log(4*4)+log(3*4))/4 rtol=1e-5
    end

    @testset "High P0, number of decisions five and initial influx high" begin
        P0 = 0.5
        I = 0.15
        decision_step = 5.0
        number_decision = 6
        deterministic = true
        max_options = 10

        entropy = PathwayDiversity._entropy(P0, I, decision_step, number_decision;
                                            deterministic=deterministic, max_options=max_options)

        @test entropy ≈ 11.856449 rtol=1e-5
    end
end

@testset "#_entropy closer_more_likely" begin
    @testset "Low P0, low influx and number of decisions one" begin
        P0 = 0.5
        I = 0.1
        decision_step = 5.0
        number_decision = 1
        deterministic = true
        max_options = 10
        method = "closer_more_likely"

        entropy = PathwayDiversity._entropy(P0, I, decision_step, number_decision;
                                            deterministic=deterministic, max_options=max_options, method=method)

        @test entropy ≈ (2/30)*log(30)+0.2*log(10)+0.4*log(5)+(1/3)*log(3) rtol=1e-5
    end

    @testset "Low P0, high influx and number of decisions one" begin
        P0 = 0.5
        I = 0.38
        decision_step = 5.0
        number_decision = 1
        deterministic = true
        max_options = 10
        method = "closer_more_likely"

        entropy = PathwayDiversity._entropy(P0, I, decision_step, number_decision;
                                            deterministic=deterministic, max_options=max_options, method=method)

        @test entropy ≈ (1/20)*log(20)+(3/20)*log(20/3)+(6/20)*log(10/3)+(1/2)*log(2) rtol=1e-5
    end

    @testset "High P0, low influx and number of decisions one" begin
        P0 = 2.0
        I = 0.06
        decision_step = 5.0
        number_decision = 1
        deterministic = true
        max_options = 10
        method = "closer_more_likely"

        entropy = PathwayDiversity._entropy(P0, I, decision_step, number_decision;
                                            deterministic=deterministic, max_options=max_options, method=method)

        @test entropy ≈ (3/25)*log(25/3)+(2*6/25)*log(25/6)+(10/25)*log(25/10) rtol=1e-5
    end

    @testset "High P0, high influx and number of decisions one" begin
        P0 = 2.0
        I = 0.38
        decision_step = 5.0
        number_decision = 1
        deterministic = true
        max_options = 10
        method = "closer_more_likely"

        entropy = PathwayDiversity._entropy(P0, I, decision_step, number_decision;
                                            deterministic=deterministic, max_options=max_options, method=method)
        @test entropy ≈ (1/20)*log(20)+(3/20)*log(20/3)+(3/10)*log(10/3)+(1/2)*log(2) rtol=1e-5
    end

    # TODO: Add tests of max_options, minimum_influx and maximun_influx
end
