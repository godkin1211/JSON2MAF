"""
runtests.jl

JSON2MAF Test Suite Main Entry Point
"""

using Test
using JSON2MAF

@testset "JSON2MAF Test Suite" begin
    @testset "Data Structure Tests" begin
        include("datastructures_tests.jl")
    end

    @testset "JSON Parsing Tests" begin
        include("parser_tests.jl")
    end

    @testset "Quality Filter Tests" begin
        include("quality_filter_tests.jl")
    end

    @testset "ClinVar Filter Tests" begin
        include("clinvar_filter_tests.jl")
    end

    @testset "Predictive Scores Tests" begin
        include("predictive_scores_tests.jl")
    end

    @testset "Integrated Decision Engine Tests" begin
        include("decision_engine_tests.jl")
    end

    @testset "MAF Conversion Tests" begin
        include("maf_converter_tests.jl")
    end

    @testset "MAF Writer Tests" begin
        include("maf_writer_tests.jl")
    end
end
