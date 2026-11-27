"""Asynchronous, cached batch access to the KEGG REST API.

The functions here mirror their synchronous counterparts in
:mod:`kegger.kegg_tools`, but fetch many records concurrently behind a
persistent HTTP cache. Every public function comes in two flavours:

* a blocking version (``org_annotation``, ``get_paths``) for scripts, and
* an ``_async`` version for callers that already run an event loop,
  such as Jupyter notebooks.
"""

import asyncio
import io
import logging
from datetime import timedelta
from typing import Any

import pandas as pd
from aiohttp import ClientTimeout
from aiohttp_client_cache import CachedSession, SQLiteBackend

logger = logging.getLogger(__name__)

#: KEGG asks clients to stay around three requests per second.
DEFAULT_MAX_CONCURRENT = 3

#: Separate from the ``requests-cache`` database used by :mod:`kegger.kegg_tools`;
#: the two libraries do not share a schema.
DEFAULT_CACHE_NAME = "kegg_async_cache.sqlite"

DEFAULT_EXPIRE_AFTER = timedelta(days=90)

#: Column names shared with :func:`kegger.kegg_tools.get_org`.
ORG_COLUMNS = ["gene", "feature", "position", "annotation"]


async def fetch_url(
    session: CachedSession,
    url: str,
    semaphore: asyncio.Semaphore,
    timeout: float = 15,
) -> dict[str, Any]:
    """Fetch a single KEGG URL through the shared cached session.

    Args:
        session: The cached session used for the request.
        url: Fully formatted KEGG REST URL.
        semaphore: Limits how many requests are in use at once.
        timeout: Total per-request timeout in seconds.

    Returns:
        dict: Always the same shape, so callers never have to guess which
        keys are present:

        - 'url': The requested URL.
        - 'raw': The response body, or None if the request failed.
        - 'error': None on success, otherwise a short description.
        - 'from_cache': Whether the response was served from the local cache.
        - 'created_on': When a cached response was stored (None if live).
        - 'expires_on': When a cached response expires (None if live).
    """
    async with semaphore:
        async with session.get(url, timeout=ClientTimeout(total=timeout)) as response:
            result = {
                "url": url,
                "raw": None,
                "error": None,
                "from_cache": response.from_cache,
                "created_on": response.created_at,
                "expires_on": response.expires,
            }
            if response.status == 200:
                result["raw"] = await response.text(encoding="utf-8")
            elif response.status == 404:
                result["error"] = "HTTP 404 Not Found"
                logger.warning("Url for: %s, 404 error Not Found", url)
            else:
                result["error"] = f"HTTP {response.status}"
                logger.warning("Url for: %s, %s", url, result["error"])
            return result


async def fetch_all(
    query: list[str],
    url: str,
    max_concurrent: int = DEFAULT_MAX_CONCURRENT,
    cache_name: str = DEFAULT_CACHE_NAME,
    expire_after: timedelta = DEFAULT_EXPIRE_AFTER,
) -> list[dict[str, Any]]:
    """Fetch one KEGG URL per query item, concurrently and in order.

    Args:
        query: Items to substitute into ``url``, one request each.
        url: URL template containing a single ``{}`` placeholder.
        max_concurrent: Maximum simultaneous requests. Keep this low; KEGG
            returns HTTP 403 to clients that hammer the API.
        cache_name: Path to the SQLite cache database.
        expire_after: How long cached responses stay valid.

    Returns:
        list[dict]: One result per query item, in the same order as ``query``.
        Requests that raised (timeouts, connection errors) come back as
        error dicts rather than exceptions.
    """
    semaphore = asyncio.Semaphore(max_concurrent)
    cache = SQLiteBackend(cache_name=cache_name, expire_after=expire_after)

    async with CachedSession(cache=cache) as session:
        tasks = [
            asyncio.create_task(fetch_url(session, url.format(q), semaphore))
            for q in query
        ]
        raw_results = await asyncio.gather(*tasks, return_exceptions=True)

    results = []
    for q, result in zip(query, raw_results):
        if isinstance(result, BaseException):
            logger.warning("Request for %s failed: %r", q, result)
            results.append(
                {
                    "url": url.format(q),
                    "raw": None,
                    "error": f"{type(result).__name__}: {result}",
                    "from_cache": False,
                    "created_on": None,
                    "expires_on": None,
                }
            )
        else:
            results.append(result)
    return results


def _run_sync(coro):
    """Run a coroutine, with a useful error when a loop is already running."""
    try:
        asyncio.get_running_loop()
    except RuntimeError:
        return asyncio.run(coro)

    coro.close()
    raise RuntimeError(
        "An event loop is already running (e.g. inside Jupyter). "
        "Await the async version of this function instead: "
        "org_annotation_async() or get_paths_async()."
    )


def kegg_to_dataframe(response: str, col_names: list[str]) -> pd.DataFrame:
    """Parse a tab-delimited KEGG response into a DataFrame.

    Args:
        response: Raw response body from a KEGG REST call.
        col_names: Column names to apply to the headerless table.

    Returns:
        pd.DataFrame: The parsed table, empty if the response had no rows.
    """
    if not response or not response.strip():
        return pd.DataFrame(columns=col_names)
    record = io.StringIO(response)
    return pd.read_csv(record, sep="\t", header=None, names=col_names)


def _verbose(result: dict[str, Any], func_name: str) -> None:
    """Print cache metadata for a single result."""
    print(f"{func_name}")
    print("\tUrl:", result.get("url"))
    print("\tFrom cache:", result.get("from_cache"))
    print("\tCreated on:", result.get("created_on"))
    print("\tExpires on:", result.get("expires_on"))


async def org_annotation_async(
    org_list: list[str],
    verbose: bool = False,
    max_concurrent: int = DEFAULT_MAX_CONCURRENT,
) -> pd.DataFrame:
    """Retrieve gene-level annotation for several organisms at once.

    The async counterpart of :func:`org_annotation`; see it for details.
    """
    url = "https://rest.kegg.jp/list/{}"
    results = await fetch_all(org_list, url, max_concurrent=max_concurrent)

    org_info = []
    for org, result in zip(org_list, results):
        if result["error"] is not None:
            continue
        if verbose:
            print(f"***** Org Info: {org}")
            _verbose(result, "org_annotation")

        df = kegg_to_dataframe(result["raw"], ORG_COLUMNS)
        if df.empty:
            logger.warning("No annotation rows returned for %s", org)
            continue
        # Removing the 'org' prefix, matching kegg_tools.get_org
        df["gene"] = df["gene"].str.replace(f"{org}:", "", regex=False)
        df["org"] = org
        org_info.append(df)

    if not org_info:
        logger.warning("No annotations retrieved for any of: %s", ", ".join(org_list))
        return pd.DataFrame(columns=ORG_COLUMNS + ["org"])
    return pd.concat(org_info, ignore_index=True)


def org_annotation(
    org_list: list[str],
    verbose: bool = False,
    max_concurrent: int = DEFAULT_MAX_CONCURRENT,
) -> pd.DataFrame:
    """Retrieve gene-level annotation for several organisms at once.

    Args:
        org_list: KEGG organism codes, e.g. ``["eco", "shn"]``.
        verbose: Print cache metadata for each organism.
        max_concurrent: Maximum simultaneous requests.

    Returns:
        pd.DataFrame: The concatenated annotation tables, with the same
        columns as :func:`kegger.kegg_tools.get_org` plus an 'org' column
        identifying the source organism. Organisms that failed to download
        are skipped and logged; if none succeed, the frame is empty.

    Raises:
        RuntimeError: If called while an event loop is running. Use
            :func:`org_annotation_async` from a notebook instead.

    Example:
        >>> from kegger.batch_tools import org_annotation
        >>> df = org_annotation(["eco", "shn"])
    """
    return _run_sync(org_annotation_async(org_list, verbose, max_concurrent))


async def get_paths_async(
    pathids: list[str],
    verbose: bool = False,
    max_concurrent: int = DEFAULT_MAX_CONCURRENT,
) -> list[dict[str, Any]]:
    """Download several KEGG records at once.

    The async counterpart of :func:`get_paths`; see it for details.
    """
    url = "https://rest.kegg.jp/get/{}"
    results = await fetch_all(pathids, url, max_concurrent=max_concurrent)

    results_list = []
    for pathid, result in zip(pathids, results):
        if result["error"] is not None:
            logger.warning("Skipping %s: %s", pathid, result["error"])
            continue
        if verbose:
            _verbose(result, "get_paths")
        results_list.append(result)
    return results_list


def get_paths(
    pathids: list[str],
    verbose: bool = False,
    max_concurrent: int = DEFAULT_MAX_CONCURRENT,
) -> list[dict[str, Any]]:
    """Download several KEGG records at once.

    Args:
        pathids: KEGG entry identifiers, e.g. ``["shn03020", "shn02024"]``.
        verbose: Print cache metadata for each record.
        max_concurrent: Maximum simultaneous requests.

    Returns:
        list[dict]: One result dict per successfully downloaded entry, in
        request order. The flat-file text is under the 'raw' key and can be
        passed to :func:`kegger.kegg_tools.kegg_parser`. Failed requests are
        skipped and logged, so this list may be shorter than ``pathids``.

    Raises:
        RuntimeError: If called while an event loop is running. Use
            :func:`get_paths_async` from a notebook instead.

    Example:
        >>> from kegger.batch_tools import get_paths
        >>> from kegger import kegg_parser
        >>> records = [kegg_parser(r["raw"]) for r in get_paths(["shn03020"])]
    """
    return _run_sync(get_paths_async(pathids, verbose, max_concurrent))
