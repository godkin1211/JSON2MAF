"""
clinvar_filter_tests.jl

Test ClinVar pathogenicity assessment and priority judgment functionality
"""

using Test
using JSON2MAF

@testset "ClinVar Filter Tests" begin

    @testset "Empty Entry Handling" begin
        # Test empty ClinVar entry list
        assessment = assess_clinvar_pathogenicity(ClinVarEntry[])
        @test assessment.is_pathogenic == false
        @test assessment.is_likely_pathogenic == false
        @test assessment.selected_entry === nothing
        @test assessment.confidence_level == "none"
        @test occursin("No ClinVar entries", assessment.reason)
    end

    @testset "Single Pathogenic Annotation" begin
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

    @testset "Single Likely Pathogenic Annotation" begin
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

    @testset "Benign/Uncertain Significance Filtering" begin
        # Benign entries should be filtered
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

        # Uncertain significance entries should be filtered
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

    @testset "Review Status Priority Sorting" begin
        # Test get_review_status_priority function
        @test get_review_status_priority("practice guideline") == 1
        @test get_review_status_priority("reviewed by expert panel") == 2
        @test get_review_status_priority("criteria provided, multiple submitters, no conflicts") == 3
        @test get_review_status_priority("criteria provided, conflicting interpretations") == 4
        @test get_review_status_priority("criteria provided, single submitter") == 5
        @test get_review_status_priority("no assertion criteria provided") == 6
        @test get_review_status_priority("no assertion provided") == 7
        @test get_review_status_priority("unknown status") == 8
    end

    @testset "Multiple Annotations - Review Status Priority" begin
        # Low priority entry (single submitter)
        entry_low = ClinVarEntry(
            "RCV11111", "11111",
            "Pathogenic",
            "criteria provided, single submitter",
            ["Disease A"],
            "2024-03-01"
        )

        # High priority entry (expert panel)
        entry_high = ClinVarEntry(
            "RCV22222", "22222",
            "Likely pathogenic",
            "reviewed by expert panel",
            ["Disease B"],
            "2024-01-01"
        )

        assessment = assess_clinvar_pathogenicity([entry_low, entry_high])

        # Should select the expert panel entry
        @test assessment.selected_entry === entry_high
        @test assessment.confidence_level == "high"
        @test occursin("Selected from 2 entries", assessment.reason)
    end

    @testset "Multiple Annotations - Cancer-Related Priority" begin
        # Same review status, but one is cancer-related
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

        # Should select the cancer-related entry
        @test assessment.selected_entry === entry_cancer
        @test occursin("Cancer-related disease", assessment.reason)
    end

    @testset "Multiple Annotations - Most Recent Time Priority" begin
        # Same review status and disease type, but different update times
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

        # Should select the newer entry
        @test assessment.selected_entry === entry_new
    end

    @testset "Cancer Keyword Detection" begin
        # Test various cancer keywords
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

        # Test non-cancer diseases
        @test is_cancer_related(["Diabetes"]) == false
        @test is_cancer_related(["Hypertension"]) == false
        @test is_cancer_related(["Alzheimer disease"]) == false

        # Test case insensitivity
        @test is_cancer_related(["BREAST CANCER"]) == true
        @test is_cancer_related(["Lung Carcinoma"]) == true

        # Test multiple disease names
        @test is_cancer_related(["Diabetes", "Lung cancer", "Hypertension"]) == true
        @test is_cancer_related(["Diabetes", "Hypertension"]) == false
    end

    @testset "Complex Scenario - Multiple Priority Combinations" begin
        # Scenario: 3 entries with different review status, cancer relevance, and time
        entry1 = ClinVarEntry(
            "RCV11111", "11111",
            "Pathogenic",
            "practice guideline",
            ["Diabetes"],  # Non-cancer
            "2023-01-01"   # Old time
        )

        entry2 = ClinVarEntry(
            "RCV22222", "22222",
            "Likely pathogenic",
            "criteria provided, single submitter",  # Low priority
            ["Lung cancer"],  # Cancer-related
            "2024-01-01"      # New time
        )

        entry3 = ClinVarEntry(
            "RCV33333", "33333",
            "Pathogenic",
            "reviewed by expert panel",  # Medium priority
            ["Breast carcinoma"],  # Cancer-related
            "2024-01-01"           # New time
        )

        assessment = assess_clinvar_pathogenicity([entry1, entry2, entry3])

        # Should select entry1 because practice guideline has highest priority
        @test assessment.selected_entry === entry1
        @test assessment.confidence_level == "high"
    end

    @testset "Confidence Level Assessment" begin
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

    @testset "Missing Field Handling" begin
        # Test when last_evaluated is nothing
        entry1 = ClinVarEntry(
            "RCV11111", "11111",
            "Pathogenic",
            "criteria provided, single submitter",
            ["Disease"],
            nothing  # Missing time
        )

        entry2 = ClinVarEntry(
            "RCV22222", "22222",
            "Pathogenic",
            "criteria provided, single submitter",
            ["Disease"],
            "2024-01-01"
        )

        # Should handle normally, prioritize one with time
        assessment = assess_clinvar_pathogenicity([entry1, entry2])
        @test assessment.selected_entry === entry2
    end

    @testset "Pathogenic/Likely Pathogenic Mixed Assessment" begin
        # Test "Pathogenic/Likely pathogenic" combination
        entry = ClinVarEntry(
            "RCV11111", "11111",
            "Pathogenic/Likely pathogenic",
            "criteria provided, multiple submitters, no conflicts",
            ["Disease"],
            "2024-01-01"
        )

        # This should be identified as pathogenic (contains pathogenic but not just likely)
        assessment = assess_clinvar_pathogenicity([entry])
        @test assessment.is_pathogenic == true
        @test assessment.is_likely_pathogenic == false
    end
end
