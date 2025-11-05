"""
runtests.jl

JSON2MAF 測試套件主入口
"""

using Test
using JSON2MAF

@testset "JSON2MAF 測試套件" begin
    @testset "資料結構測試" begin
        include("datastructures_tests.jl")
    end

    @testset "JSON 解析測試" begin
        include("parser_tests.jl")
    end

    @testset "品質過濾測試" begin
        include("quality_filter_tests.jl")
    end

    @testset "ClinVar 過濾測試" begin
        include("clinvar_filter_tests.jl")
    end

    @testset "預測分數測試" begin
        include("predictive_scores_tests.jl")
    end

    @testset "整合決策引擎測試" begin
        include("decision_engine_tests.jl")
    end

    @testset "MAF 轉換測試" begin
        include("maf_converter_tests.jl")
    end

    @testset "MAF 寫入測試" begin
        include("maf_writer_tests.jl")
    end
end
