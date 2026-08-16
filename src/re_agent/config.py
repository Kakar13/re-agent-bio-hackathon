from pydantic_settings import BaseSettings, SettingsConfigDict


class Settings(BaseSettings):
    model_config = SettingsConfigDict(
        env_file=".env",
        env_file_encoding="utf-8",
        extra="ignore",
    )

    anthropic_api_key: str = ""
    proto_api_key: str = ""
    hf_token: str = ""
    modal_token_id: str = ""
    modal_token_secret: str = ""

    langsmith_api_key: str = ""
    langsmith_tracing: str = ""
    langsmith_project: str = "reAgent-hackathon"

    benchling_tenant_url: str = ""
    benchling_url: str = ""
    benchling_api_key: str = ""
    benchling_folder_id: str = ""
    benchling_run_schema_id: str = ""

    netmhcpan_bin: str = ""
    netmhciipan_bin: str = ""

    def missing(self) -> list[str]:
        checks = {
            "ANTHROPIC_API_KEY": self.anthropic_api_key,
            "PROTO_API_KEY": self.proto_api_key,
        }
        return [name for name, value in checks.items() if not value]


settings = Settings()
