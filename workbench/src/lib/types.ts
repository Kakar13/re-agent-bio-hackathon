export type Message = {
  id?: string;
  type?: string;
  role?: string;
  tool_call_id?: string;
  name?: string;
  content: string | Array<{ type: string; text?: string }>;
  tool_calls?: Array<{
    id: string;
    name: string;
    args: Record<string, unknown>;
  }>;
};

export type ReviewCheck = {
  name: string;
  passed: boolean;
  detail: string;
};

export type Review = {
  artifact_id: string;
  status: "pass" | "fail";
  checks: ReviewCheck[];
  summary: string;
};

export type Assessment = {
  candidate_id: string;
  sequence: string;
  combined_rank_score: number | null;
  response_results: Array<{
    adapter_id: string;
    status: string;
    score_scale: string;
  }>;
  mhc_i_surrogate_results: Array<{
    adapter_id: string;
    status: string;
    allele: string;
    protein_summary: Record<string, number>;
    spatial_tracks: Record<string, number[]>;
    provenance: {
      provider: string;
      version: string;
      capability: string;
      source: string;
      parameters: {
        checkpoint_sha256?: string;
        dataset_version_hash?: string;
        encoder_model_id?: string;
        [key: string]: unknown;
      };
      input_sha256: string;
      runtime_seconds: number;
      cached: boolean;
      timestamp: string;
    };
    predictions: Array<{
      start: number;
      end: number;
      peptide: string;
      cleavage_n_probability: number;
      cleavage_c_probability: number;
      tap_log_ic50_relative: number;
      tap_uncertainty: number;
      mhc_i_presentation_propensity: number;
      composite_processing_risk: number;
      confidence: number;
    }>;
  }>;
  mhc_results: Array<{
    provider_id: string;
    status: string;
    supported_alleles: string[];
    hits: Array<{
      allele: string;
      start: number;
      end: number;
      peptide: string;
      core: string;
      el_rank?: number;
      ba_rank?: number;
      ba_ic50_nm?: number;
    }>;
    provenance: Record<string, unknown>;
  }>;
  processing: Array<Record<string, unknown>>;
  tolerance: Array<Record<string, unknown>>;
  spatial_tracks: Record<string, number[]>;
  structure?: {
    path: string;
    format: "pdb";
    chain_id: string;
    residue_ids: string[];
    unresolved_sequence_positions: number[];
    sequence_sha256: string;
    structure_sha256: string;
    mapping_status: "verified_exact_sequence" | "verified_terminal_trim";
  };
  warnings: string[];
};

export type ScientificArtifact = {
  id: string;
  kind: string;
  title: string;
  path: string;
  manifest_path: string;
  sha256: string;
  payload: {
    title: string;
    claim_boundary?: string;
    citations?: Array<{ claim: string; url: string }>;
    assessment?: Assessment;
    architecture?: {
      status: string;
      checks: Record<string, { available: boolean }>;
      providers: {
        netmhciipan: { channels: string[] };
        response_model: { available: boolean };
      };
    };
    manifest?: {
      status: string;
      spec_path: string;
      candidate_counts: {
        generated?: number;
        validation_passed?: number;
        validation_failed?: number;
        screened?: number;
      };
      candidates: Array<{
        candidate_id: string;
        validation_status: string;
        screening_status: string;
        validation_checks: Array<{
          name: string;
          recorded: boolean;
          passed: boolean;
          value?: string | number | boolean | null;
        }>;
      }>;
    };
    assessments?: Assessment[];
    spec_path?: string;
    status?: string;
    checks?: ReviewCheck[];
    summary?: string;
  };
};

export type AgentState = {
  [key: string]: unknown;
  messages: Message[];
  artifacts: ScientificArtifact[];
  reviews: Review[];
  screening_profile?:
    | "mhc_ii_standard"
    | "mhc_ii_plus_chao1";
  direct_screen_request?: {
    sequence: string;
    candidate_id: string;
  } | null;
};
