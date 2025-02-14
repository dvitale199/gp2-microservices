from pathlib import Path
from pydantic_settings import BaseSettings

class Settings(BaseSettings):
    PROJECT_ROOT: Path = Path(__file__).parent.parent
    PLINK2_PATH: str = str(PROJECT_ROOT / "exec" / "plink2")
    PROJECT_DIR: str = '/home/vitaled2/gp2_carriers'
    GCLOUD_BUCKET: str = 'gs://gp2_carriers/'
    
    class Config:
        env_file = ".env"

settings = Settings()