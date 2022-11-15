import importlib.metadata
import os
from datetime import datetime, timezone
from pathlib import Path
from typing import *

import mechanicalsoup
from platformdirs import PlatformDirs
from typer_cloup import *


_T = TypeVar("_T")

GenericPath = Union[AnyStr, os.PathLike[AnyStr]]

package = importlib.metadata.metadata(__package__)
platform_dirs = PlatformDirs(appname = package["Name"], version = package["Version"])


def first(seq: Sequence[_T]) -> Optional[_T]:
	return seq[0] if len(seq) > 0 else None


def utc_to_local(utc_dt: datetime) -> datetime:
	return utc_dt.replace(tzinfo = timezone.utc).astimezone(tz = None)


def debug_mode() -> bool:
	return os.getenv("DEBUG") == "1"


def cache_dir() -> Path:
	path = Path(platform_dirs.user_cache_dir)
	path.mkdir(parents = True, exist_ok = True)
	return path


def _form_choose_submit_none(form: mechanicalsoup.Form):
	buttons = (button for button in form.form.select("input[type='submit' i], button") if button.get("type", "").lower() not in ("button", "reset"))

	for button in buttons:
		del button["name"]

	form._submit_chosen = True


mechanicalsoup.Form.choose_submit_none = _form_choose_submit_none
