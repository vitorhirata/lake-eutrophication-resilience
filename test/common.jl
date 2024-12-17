@testset "#entropy standard calculus" begin
    @testset "Low P0 and number of decisions one" begin
        P0 = 0.5
        I = 0.1
        decision_step = 5.0
        number_decision = 1
        max_options = 10

        entropy = PathwayDiversity.entropy(P0, I, decision_step, number_decision)
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
        max_options = 10

        entropy = PathwayDiversity.entropy(P0, I, decision_step, number_decision)
        number_options = PathwayDiversity.number_possible_influx(P0, max_options)

        @test number_options == 4
        @test entropy ≈ log(number_options) rtol=1e-5
    end

    @testset "Low P0 and number of decisions two" begin
        P0 = 0.5
        I = 0.1
        decision_step = 5.0
        number_decision = 2

        entropy = PathwayDiversity.entropy(P0, I, decision_step, number_decision)

        # To compute the entropy in this case the number of influx from the last state will simplify with the
        # probability of choosing the influx. So, from each final state we will have:
        # last_n_opt*(1/first_n_opt*last_n_opt)*log(1/first_n_opt*last_n_opt)
        # Ex. -9*(1/9*10)*log(1/90) = (1/10)*log(90)
        @test entropy ≈ (5*log(100)+2*log(90)+2*log(80)+log(70))/10 rtol=1e-5
    end

    @testset "High P0 and number of decisions two" begin
        P0 = 2.0
        I = 0.1
        decision_step = 5.0
        number_decision = 2

        entropy = PathwayDiversity.entropy(P0, I, decision_step, number_decision)

        @test entropy ≈ (2*log(4*4)+2*log(3*4))/4 rtol=1e-5
    end

    @testset "High P0, number of decisions five and initial influx high" begin
        P0 = 0.5
        I = 0.15
        decision_step = 5.0
        number_decision = 6

        entropy = PathwayDiversity.entropy(P0, I, decision_step, number_decision)

        @test entropy ≈ 12.492185 rtol=1e-5
    end

    @testset "Lower P0, number of decisions one and max_options 15" begin
        P0 = 0.2
        I = 0.1
        decision_step = 5.0
        number_decision = 1
        max_options = 15

        entropy = PathwayDiversity.entropy(P0, I, decision_step, number_decision; max_options=max_options)
        number_options = PathwayDiversity.number_possible_influx(P0, max_options)

        @test number_options == 15
        @test entropy ≈ log(number_options) rtol=1e-5
    end

    @testset "Low P0, number of decisions two and minimum_influx lower" begin
        P0 = 0.5
        I = 0.1
        decision_step = 5.0
        number_decision = 2
        minimum_influx = 0.02

        entropy = PathwayDiversity.entropy(P0, I, decision_step, number_decision; minimum_influx=minimum_influx)

        @test entropy ≈ (6*log(100)+log(90)+2*log(80)+log(70))/10 rtol=1e-5
    end

    @testset "Low P0, number of decisions two and maximum_influx higher" begin
        P0 = 0.5
        I = 0.1
        decision_step = 5.0
        number_decision = 2
        maximum_influx = 0.32

        entropy = PathwayDiversity.entropy(P0, I, decision_step, number_decision; maximum_influx=maximum_influx)

        @test entropy ≈ (5*log(100)+2*log(90)+log(80)+log(70)+log(60))/10 rtol=1e-5
    end
end

@testset "#entropy deterministic false" begin
    @testset "Low P0 and number of decisions two" begin
        P0 = 0.5
        I = 0.1
        decision_step = 5.0
        number_decision = 2
        deterministic = false

        Random.seed!(1234)

        entropy = PathwayDiversity.entropy(P0, I, decision_step, number_decision; deterministic=deterministic)

        @test entropy ≈ (3*log(90)+2*log(80)+log(70)+log(60)+2*log(50)+log(40))/10 rtol=1e-5
    end

    @testset "High P0 and number of decisions two" begin
        P0 = 2.0
        I = 0.1
        decision_step = 5.0
        number_decision = 2
        deterministic = false

        Random.seed!(1234)

        entropy = PathwayDiversity.entropy(P0, I, decision_step, number_decision; deterministic=deterministic)

        @test entropy ≈ (3*log(2*4)+log(1*4))/4 rtol=1e-5
    end
end

@testset "#entropy closer_more_likely" begin
    @testset "Low P0, low influx and number of decisions one" begin
        P0 = 0.5
        I = 0.1
        decision_step = 5.0
        number_decision = 1
        method = "closer_more_likely"

        entropy = PathwayDiversity.entropy(P0, I, decision_step, number_decision; method=method)

        @test entropy ≈ (2/30)*log(30)+0.2*log(10)+0.4*log(5)+(1/3)*log(3) rtol=1e-5
    end

    @testset "Low P0, high influx and number of decisions one" begin
        P0 = 0.5
        I = 0.38
        decision_step = 5.0
        number_decision = 1
        method = "closer_more_likely"

        entropy = PathwayDiversity.entropy(P0, I, decision_step, number_decision; method=method)

        @test entropy ≈ (1/20)*log(20)+(3/20)*log(20/3)+(6/20)*log(10/3)+(1/2)*log(2) rtol=1e-5
    end

    @testset "High P0, low influx and number of decisions one" begin
        P0 = 2.0
        I = 0.06
        decision_step = 5.0
        number_decision = 1
        method = "closer_more_likely"

        entropy = PathwayDiversity.entropy(P0, I, decision_step, number_decision; method=method)

        @test entropy ≈ (3/25)*log(25/3)+(2*6/25)*log(25/6)+(10/25)*log(25/10) rtol=1e-5
    end

    @testset "High P0, high influx and number of decisions one" begin
        P0 = 2.0
        I = 0.38
        decision_step = 5.0
        number_decision = 1
        method = "closer_more_likely"

        entropy = PathwayDiversity.entropy(P0, I, decision_step, number_decision; method=method)
        @test entropy ≈ (1/20)*log(20)+(3/20)*log(20/3)+(3/10)*log(10/3)+(1/2)*log(2) rtol=1e-5
    end
end

@testset "#entropy further_more_likely" begin
    @testset "Low P0, low influx and number of decisions one" begin
        P0 = 0.5
        I = 0.1
        decision_step = 5.0
        number_decision = 1
        method = "further_more_likely"

        entropy = PathwayDiversity.entropy(P0, I, decision_step, number_decision; method=method)

        @test entropy ≈ log(63)-(2*7*log(7)+2*4*log(4)+2*1*log(1)+10*log(10)+13*log(13)+16*log(16))/63 rtol=1e-5
    end

    @testset "Low P0, high influx and number of decisions one" begin
        P0 = 0.5
        I = 0.38
        decision_step = 5.0
        number_decision = 1
        method = "further_more_likely"

        entropy = PathwayDiversity.entropy(P0, I, decision_step, number_decision; method=method)

        @test entropy ≈ log(51)-(16*log(16)+13*log(13)+10*log(10)+7*log(7)+4*log(4)+1*log(1))/51 rtol=1e-5
    end

    @testset "High P0, low influx and number of decisions one" begin
        P0 = 2.0
        I = 0.06
        decision_step = 5.0
        number_decision = 1
        method = "further_more_likely"

        entropy = PathwayDiversity.entropy(P0, I, decision_step, number_decision; method=method)

        @test entropy ≈ log(52)-(16*log(16)+2*13*log(13)+10*log(10))/52 rtol=1e-5
    end

    @testset "High P0, high influx and number of decisions one" begin
        P0 = 2.0
        I = 0.38
        decision_step = 5.0
        number_decision = 1
        method = "further_more_likely"

        entropy = PathwayDiversity.entropy(P0, I, decision_step, number_decision; method=method)
        @test entropy ≈ log(46)-(16*log(16)+13*log(13)+10*log(10)+7*log(7))/46 rtol=1e-5
    end
end
