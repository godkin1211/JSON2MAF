"""
parser_tests.jl

測試 Nirvana JSON 解析功能
"""

using Test
using JSON2MAF

@testset "Nirvana Parser 測試" begin
    # 使用相對於專案根目錄的路徑
    project_root = dirname(@__DIR__)
    test_file = joinpath(project_root, "examples", "sample_data", "test_mini.json.gz")

    @testset "檔案讀取測試" begin
        @test isfile(test_file)

        # 測試不存在的檔案
        @test_throws ErrorException parse_nirvana_json("nonexistent.json.gz")
    end

    @testset "完整解析測試" begin
        data = parse_nirvana_json(test_file)

        @test data isa NirvanaData
        @test data.header isa NirvanaHeader
        @test length(data.positions) > 0
    end

    @testset "Header 解析測試" begin
        data = parse_nirvana_json(test_file)
        header = data.header

        @test header.annotator == "Test Annotator 1.0"
        @test header.creation_time == "2025-10-30 12:00:00"
        @test header.genome_assembly == "GRCh38"
        @test header.schema_version == 6
        @test length(header.data_sources) == 1
        @test header.data_sources[1].name == "ClinVar"
        @test length(header.samples) == 1
        @test header.samples[1] == "TEST_SAMPLE"
    end

    @testset "Position 解析測試" begin
        data = parse_nirvana_json(test_file)
        @test length(data.positions) == 1

        pos = data.positions[1]

        # 基本資訊
        @test pos.chromosome == "chr1"
        @test pos.start == 12345
        @test pos.reference_allele == "A"
        @test pos.alternate_allele == "G"
        @test pos.variant_type == "SNV"

        # 樣本資訊
        @test pos.total_depth == 100
        @test pos.variant_frequencies !== nothing
        @test length(pos.variant_frequencies) == 1
        @test pos.variant_frequencies[1] == 0.45
    end

    @testset "ClinVar 解析測試" begin
        data = parse_nirvana_json(test_file)
        pos = data.positions[1]

        @test length(pos.clinvar) == 1
        cv = pos.clinvar[1]

        @test cv.id == "RCV000001"
        @test cv.clinical_significance == "Pathogenic"
        @test cv.review_status == "criteria provided, multiple submitters"
        @test length(cv.diseases) == 1
        @test cv.diseases[1] == "Cancer"
        @test cv.last_evaluated == "2024-01-01"
    end

    @testset "COSMIC 解析測試" begin
        data = parse_nirvana_json(test_file)
        pos = data.positions[1]

        @test length(pos.cosmic) >= 1
        cosmic = pos.cosmic[1]

        @test cosmic.id == "COSV12345"
        @test cosmic.count !== nothing
    end

    @testset "Transcript 解析測試" begin
        data = parse_nirvana_json(test_file)
        pos = data.positions[1]

        @test length(pos.transcripts) == 1
        trans = pos.transcripts[1]

        @test trans.id == "NM_001234.1"
        @test trans.gene_symbol == "BRCA1"
        @test length(trans.consequence) == 1
        @test trans.consequence[1] == "missense_variant"
        @test trans.hgvsc == "c.123A>G"
        @test trans.hgvsp == "p.Thr41Ala"
    end

    @testset "族群頻率解析測試" begin
        data = parse_nirvana_json(test_file)
        pos = data.positions[1]

        @test length(pos.population_frequencies) >= 2

        # 查找 gnomad-exome
        gnomad_exome = findfirst(pf -> pf.source == "gnomad-exome", pos.population_frequencies)
        @test gnomad_exome !== nothing

        gne = pos.population_frequencies[gnomad_exome]
        @test gne.all_af == 0.001
        @test gne.eas_af == 0.002

        # 查找 oneKg
        onekg = findfirst(pf -> pf.source == "oneKg", pos.population_frequencies)
        @test onekg !== nothing

        okg = pos.population_frequencies[onekg]
        @test okg.all_af == 0.003
        @test okg.eas_af == 0.004
    end

    @testset "預測分數解析測試" begin
        data = parse_nirvana_json(test_file)
        pos = data.positions[1]

        # DANN score
        @test pos.dann_score !== nothing
        @test pos.dann_score == 0.98

        # REVEL score
        @test pos.revel_score !== nothing
        @test pos.revel_score == 0.85

        # PrimateAI-3D score
        @test pos.primate_ai_3d !== nothing
        @test pos.primate_ai_3d == 0.92
    end

    @testset "dbSNP 解析測試" begin
        data = parse_nirvana_json(test_file)
        pos = data.positions[1]

        @test length(pos.dbsnp_ids) == 1
        @test pos.dbsnp_ids[1] == "rs12345"
    end

    @testset "缺失欄位處理測試" begin
        # 這個測試確保當某些欄位缺失時不會crash
        # 使用實際檔案進行測試會更好，但先用minimal檔案確保基本功能
        data = parse_nirvana_json(test_file)
        @test data isa NirvanaData

        # 測試 genes 欄位缺失
        @test data.genes isa Dict
    end
end

# 如果有實際的大型檔案，可以進行額外測試
# 由於實際檔案很大（>40MB），默認跳過此測試
# 設定環境變數 TEST_REAL_FILES=1 來啟用此測試
@testset "實際檔案解析測試 (可選)" begin
    project_root = dirname(@__DIR__)
    real_file = joinpath(project_root, "P.hard-filtered.vcf.annotated.json.gz")

    test_real_files = get(ENV, "TEST_REAL_FILES", "0") == "1"

    if isfile(real_file) && test_real_files
        @testset "實際檔案基本解析" begin
            println("  ⏳ 解析大型實際檔案 (這可能需要幾分鐘)...")
            data = parse_nirvana_json(real_file)

            @test data isa NirvanaData
            @test data.header isa NirvanaHeader
            @test length(data.positions) > 0

            println("  ℹ 實際檔案包含 $(length(data.positions)) 個變異位點")

            # 檢查第一個位點
            if length(data.positions) > 0
                pos = data.positions[1]
                @test pos.chromosome !== ""
                @test pos.start > 0

                # 顯示一些統計資訊
                has_clinvar = count(p -> length(p.clinvar) > 0, data.positions)
                has_cosmic = count(p -> length(p.cosmic) > 0, data.positions)
                println("  ℹ 有 ClinVar 註解: $has_clinvar")
                println("  ℹ 有 COSMIC 註解: $has_cosmic")
            end
        end
    else
        if !isfile(real_file)
            @test_skip "實際檔案不存在，跳過測試"
        else
            @test_skip "實際檔案測試已禁用 (設定 TEST_REAL_FILES=1 啟用)"
        end
    end
end
