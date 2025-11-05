"""
clinvar_filter_tests.jl

測試 ClinVar 致病性評估與優先級判斷功能
"""

using Test
using JSON2MAF

@testset "ClinVar 過濾器測試" begin

    @testset "空條目處理" begin
        # 測試空的 ClinVar 條目列表
        assessment = assess_clinvar_pathogenicity(ClinVarEntry[])
        @test assessment.is_pathogenic == false
        @test assessment.is_likely_pathogenic == false
        @test assessment.selected_entry === nothing
        @test assessment.confidence_level == "none"
        @test occursin("No ClinVar entries", assessment.reason)
    end

    @testset "單一 Pathogenic 註解" begin
        entry = ClinVarEntry(
            "RCV12345",         # id
            "12345",            # allele_id
            "Pathogenic",       # clinical_significance
            "criteria provided, multiple submitters, no conflicts",  # review_status
            ["Breast cancer"],  # diseases
            "2024-01-15"        # last_evaluated
        )

        assessment = assess_clinvar_pathogenicity([entry])
        @test assessment.is_pathogenic == true
        @test assessment.is_likely_pathogenic == false
        @test assessment.selected_entry === entry
        @test assessment.confidence_level == "medium"
        @test occursin("Pathogenic", assessment.reason)
        @test occursin("Cancer-related disease", assessment.reason)
    end

    @testset "單一 Likely pathogenic 註解" begin
        entry = ClinVarEntry(
            "RCV67890",         # id
            "67890",            # allele_id
            "Likely pathogenic",  # clinical_significance
            "criteria provided, single submitter",  # review_status
            ["Hereditary cancer syndrome"],  # diseases
            "2024-02-01"        # last_evaluated
        )

        assessment = assess_clinvar_pathogenicity([entry])
        @test assessment.is_pathogenic == false
        @test assessment.is_likely_pathogenic == true
        @test assessment.selected_entry === entry
        @test assessment.confidence_level == "low"
        @test occursin("Likely pathogenic", assessment.reason)
    end

    @testset "Benign/Uncertain significance 過濾" begin
        # Benign 條目應該被過濾
        benign_entry = ClinVarEntry(
            "RCV11111", "11111",
            "Benign",
            "reviewed by expert panel",
            ["Some disease"],
            "2024-01-01"
        )

        assessment = assess_clinvar_pathogenicity([benign_entry])
        @test assessment.is_pathogenic == false
        @test assessment.is_likely_pathogenic == false
        @test occursin("No pathogenic", assessment.reason)

        # Uncertain significance 條目應該被過濾
        uncertain_entry = ClinVarEntry(
            "RCV22222", "22222",
            "Uncertain significance",
            "criteria provided, single submitter",
            ["Some disease"],
            "2024-01-01"
        )

        assessment = assess_clinvar_pathogenicity([uncertain_entry])
        @test assessment.is_pathogenic == false
        @test assessment.is_likely_pathogenic == false
    end

    @testset "Review status 優先級排序" begin
        # 測試 get_review_status_priority 函數
        @test get_review_status_priority("practice guideline") == 1
        @test get_review_status_priority("reviewed by expert panel") == 2
        @test get_review_status_priority("criteria provided, multiple submitters, no conflicts") == 3
        @test get_review_status_priority("criteria provided, conflicting interpretations") == 4
        @test get_review_status_priority("criteria provided, single submitter") == 5
        @test get_review_status_priority("no assertion criteria provided") == 6
        @test get_review_status_priority("no assertion provided") == 7
        @test get_review_status_priority("unknown status") == 8
    end

    @testset "多筆註解 - Review status 優先" begin
        # 低優先級條目 (single submitter)
        entry_low = ClinVarEntry(
            "RCV11111", "11111",
            "Pathogenic",
            "criteria provided, single submitter",
            ["Disease A"],
            "2024-03-01"
        )

        # 高優先級條目 (expert panel)
        entry_high = ClinVarEntry(
            "RCV22222", "22222",
            "Likely pathogenic",
            "reviewed by expert panel",
            ["Disease B"],
            "2024-01-01"
        )

        assessment = assess_clinvar_pathogenicity([entry_low, entry_high])

        # 應該選擇 expert panel 的條目
        @test assessment.selected_entry === entry_high
        @test assessment.confidence_level == "high"
        @test occursin("Selected from 2 entries", assessment.reason)
    end

    @testset "多筆註解 - 癌症相關優先" begin
        # 相同 review status，但一個是癌症相關
        entry_non_cancer = ClinVarEntry(
            "RCV11111", "11111",
            "Pathogenic",
            "criteria provided, single submitter",
            ["Hypertension"],
            "2024-01-01"
        )

        entry_cancer = ClinVarEntry(
            "RCV22222", "22222",
            "Pathogenic",
            "criteria provided, single submitter",
            ["Colorectal carcinoma"],
            "2024-01-01"
        )

        assessment = assess_clinvar_pathogenicity([entry_non_cancer, entry_cancer])

        # 應該選擇癌症相關的條目
        @test assessment.selected_entry === entry_cancer
        @test occursin("Cancer-related disease", assessment.reason)
    end

    @testset "多筆註解 - 最新時間優先" begin
        # 相同 review status 和疾病類型，但更新時間不同
        entry_old = ClinVarEntry(
            "RCV11111", "11111",
            "Pathogenic",
            "criteria provided, single submitter",
            ["Disease A"],
            "2023-01-01"
        )

        entry_new = ClinVarEntry(
            "RCV22222", "22222",
            "Pathogenic",
            "criteria provided, single submitter",
            ["Disease A"],
            "2024-01-01"
        )

        assessment = assess_clinvar_pathogenicity([entry_old, entry_new])

        # 應該選擇較新的條目
        @test assessment.selected_entry === entry_new
    end

    @testset "癌症關鍵字檢測" begin
        # 測試各種癌症關鍵字
        @test is_cancer_related(["Breast cancer"]) == true
        @test is_cancer_related(["Lung carcinoma"]) == true
        @test is_cancer_related(["Brain tumor"]) == true
        @test is_cancer_related(["Malignant melanoma"]) == true
        @test is_cancer_related(["Colorectal neoplasm"]) == true
        @test is_cancer_related(["Non-Hodgkin lymphoma"]) == true
        @test is_cancer_related(["Acute myeloid leukemia"]) == true
        @test is_cancer_related(["Osteosarcoma"]) == true
        @test is_cancer_related(["Glioma"]) == true
        @test is_cancer_related(["Hepatoblastoma"]) == true

        # 測試非癌症疾病
        @test is_cancer_related(["Diabetes"]) == false
        @test is_cancer_related(["Hypertension"]) == false
        @test is_cancer_related(["Alzheimer disease"]) == false

        # 測試大小寫不敏感
        @test is_cancer_related(["BREAST CANCER"]) == true
        @test is_cancer_related(["Lung Carcinoma"]) == true

        # 測試多個疾病名稱
        @test is_cancer_related(["Diabetes", "Lung cancer", "Hypertension"]) == true
        @test is_cancer_related(["Diabetes", "Hypertension"]) == false
    end

    @testset "複雜情境 - 多重優先級組合" begin
        # 情境：3個條目，不同 review status、癌症相關性、時間
        entry1 = ClinVarEntry(
            "RCV11111", "11111",
            "Pathogenic",
            "practice guideline",
            ["Diabetes"],  # 非癌症
            "2023-01-01"   # 舊時間
        )

        entry2 = ClinVarEntry(
            "RCV22222", "22222",
            "Likely pathogenic",
            "criteria provided, single submitter",  # 低優先級
            ["Lung cancer"],  # 癌症相關
            "2024-01-01"      # 新時間
        )

        entry3 = ClinVarEntry(
            "RCV33333", "33333",
            "Pathogenic",
            "reviewed by expert panel",  # 中優先級
            ["Breast carcinoma"],  # 癌症相關
            "2024-01-01"           # 新時間
        )

        assessment = assess_clinvar_pathogenicity([entry1, entry2, entry3])

        # 應該選擇 entry1，因為 practice guideline 優先級最高
        @test assessment.selected_entry === entry1
        @test assessment.confidence_level == "high"
    end

    @testset "信心等級判斷" begin
        # High confidence: practice guideline, expert panel
        entry_high = ClinVarEntry(
            "RCV11111", "11111",
            "Pathogenic",
            "reviewed by expert panel",
            ["Disease"],
            "2024-01-01"
        )
        assessment = assess_clinvar_pathogenicity([entry_high])
        @test assessment.confidence_level == "high"

        # Medium confidence: multiple submitters
        entry_med = ClinVarEntry(
            "RCV22222", "22222",
            "Pathogenic",
            "criteria provided, multiple submitters, no conflicts",
            ["Disease"],
            "2024-01-01"
        )
        assessment = assess_clinvar_pathogenicity([entry_med])
        @test assessment.confidence_level == "medium"

        # Low confidence: single submitter
        entry_low = ClinVarEntry(
            "RCV33333", "33333",
            "Pathogenic",
            "criteria provided, single submitter",
            ["Disease"],
            "2024-01-01"
        )
        assessment = assess_clinvar_pathogenicity([entry_low])
        @test assessment.confidence_level == "low"
    end

    @testset "缺失欄位處理" begin
        # 測試 last_evaluated 為 nothing
        entry1 = ClinVarEntry(
            "RCV11111", "11111",
            "Pathogenic",
            "criteria provided, single submitter",
            ["Disease"],
            nothing  # 缺失時間
        )

        entry2 = ClinVarEntry(
            "RCV22222", "22222",
            "Pathogenic",
            "criteria provided, single submitter",
            ["Disease"],
            "2024-01-01"
        )

        # 應該能正常處理，有時間的優先
        assessment = assess_clinvar_pathogenicity([entry1, entry2])
        @test assessment.selected_entry === entry2
    end

    @testset "Pathogenic/Likely pathogenic 混合判斷" begin
        # 測試 "Pathogenic/Likely pathogenic" 這種組合
        entry = ClinVarEntry(
            "RCV11111", "11111",
            "Pathogenic/Likely pathogenic",
            "criteria provided, multiple submitters, no conflicts",
            ["Disease"],
            "2024-01-01"
        )

        # 這種情況下應該被識別為 pathogenic（包含 pathogenic 但不只是 likely）
        assessment = assess_clinvar_pathogenicity([entry])
        @test assessment.is_pathogenic == true
        @test assessment.is_likely_pathogenic == false
    end
end
